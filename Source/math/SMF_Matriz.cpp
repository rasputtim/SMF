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

#include "math/SMF_Matriz.h"
#include "util/SMF_StringUtils.h"
#include "math/SMF_EulerAngles.h"
#include "math/SMF_Quaternion.h"
#include "math/SMF_Rotation.h"
#include "util/SMF_Debug.h"
#include "structures/SMF_List.h"
#include "math/SMF_JointTransform.h"
#include "math/SMF_Simd.h"
//#pragma hdrstop

//TODO SEE WHERE TO PUT THIS TEMPLATES
namespace SMF{
namespace MATH{
	class CQuaternion;

template<typename Matrix>
void SetMatrixRotatePart(Matrix &m, const CQuaternion &q)
{
	// See e.g. http://www.geometrictools.com/Documentation/LinearAlgebraicQuaternions.pdf .

	SMF_ASSERT(q.isNormalized());
	const float x = q.x; const float y = q.y; const float z = q.z; const float w = q.w;
	m[0][0] = 1 - 2*(y*y + z*z); m[0][1] =     2*(x*y - z*w); m[0][2] =     2*(x*z + y*w);
	m[1][0] =     2*(x*y + z*w); m[1][1] = 1 - 2*(x*x + z*z); m[1][2] =     2*(y*z - x*w);
	m[2][0] =     2*(x*z - y*w); m[2][1] =     2*(y*z + x*w); m[2][2] = 1 - 2*(x*x + y*y);
}


/** Sets the top-left 3x3 area of the matrix to the rotation matrix about the Z-axis. Elements
	outside the top-left 3x3 area are ignored. This matrix rotates counterclockwise if multiplied
	in the order M*v, and clockwise if rotated in the order v*M.
	\param m The matrix to store the result.
	\param angle The rotation angle in radians. */
template<typename Matrix>
void Set3x3PartRotateZ(Matrix &m, float angle)
{
	/*
	 ³  cz -sz   0 ³
	 ³  sz  cz   0 ³
	 ³   0   0   1 ³
	*/

	const float cosz = CMath::cos(angle);
	const float sinz = CMath::sin(angle);

	m[0][0] = cosz; m[0][1] = -sinz; m[0][2] = 0.f;
	m[1][0] = sinz; m[1][1] =  cosz; m[1][2] = 0.f;
	m[2][0] =  0.f; m[2][1] =   0.f; m[2][2] = 1.f;
}

/** Sets the top-left 3x3 area of the matrix to the rotation matrix about the Y-axis. Elements
	outside the top-left 3x3 area are ignored. This matrix rotates counterclockwise if multiplied
	in the order M*v, and clockwise if rotated in the order v*M.
	\param m The matrix to store the result.
	\param angle The rotation angle in radians. */
template<typename Matrix>
void Set3x3PartRotateY(Matrix &m, float angle)
{
	/*
	 ³  cz   0  sz ³
	 ³   0   1   0 ³
	 ³ -sz   0  cz ³
	*/

	const float cosz = CMath::cos(angle);
	const float sinz = CMath::sin(angle);

	m[0][0] =  cosz; m[0][1] = 0.f; m[0][2] = sinz;
	m[1][0] =   0.f; m[1][1] = 1.f; m[1][2] =  0.f;
	m[2][0] = -sinz; m[2][1] = 0.f; m[2][2] = cosz;
}

/** Sets the top-left 3x3 area of the matrix to the rotation matrix about the X-axis. Elements
	outside the top-left 3x3 area are ignored. This matrix rotates counterclockwise if multiplied
	in the order M*v, and clockwise if rotated in the order v*M.
	\param m The matrix to store the result.
	\param angle The rotation angle in radians. */
template<typename Matrix>
void Set3x3PartRotateX(Matrix &m, float angle)
{
	/*
	 ³   1   0   0 ³
	 ³   0  cz -sz ³
	 ³   0  sz  cz ³
	*/

	const float cosz = CMath::cos(angle);
	const float sinz = CMath::sin(angle);

	m[0][0] = 1.f; m[0][1] =  0.f; m[0][2] =   0.f;
	m[1][0] = 0.f; m[1][1] = cosz; m[1][2] = -sinz;
	m[2][0] = 0.f; m[2][1] = sinz; m[2][2] =  cosz;
}
} //end MATH
}  //end SMF




namespace SMF {
namespace MATH{
// Static
#define MIN(x,y)     (((x) < (y)) ? (x) : (y))
#define MAX(x,y)     (((x) > (y)) ? (x) : (y))
#define MID(x,y,z)   MAX((x), MIN((y), (z)))

void PII_inverse_4x4(float* mat)
{
float d, di;
di = mat[0];
mat[0] = d = 1.0f / di;
mat[4] *= -d;
mat[8] *= -d;
mat[12] *= -d;
mat[1] *= d;
mat[2] *= d;
mat[3] *= d;
mat[5] += mat[4] * mat[1] * di;
mat[6] += mat[4] * mat[2] * di;
mat[7] += mat[4] * mat[3] * di;
mat[9] += mat[8] * mat[1] * di;
mat[10] += mat[8] * mat[2] * di;
mat[11] += mat[8] * mat[3] * di;
mat[13] += mat[12] * mat[1] * di;
mat[14] += mat[12] * mat[2] * di;
mat[15] += mat[12] * mat[3] * di;
di = mat[5];
mat[5] = d = 1.0f / di;
mat[1] *= -d;
mat[9] *= -d;
mat[13] *= -d;
mat[4] *= d;
mat[6] *= d;
mat[7] *= d;
mat[0] += mat[1] * mat[4] * di;
mat[2] += mat[1] * mat[6] * di;
mat[3] += mat[1] * mat[7] * di;
mat[8] += mat[9] * mat[4] * di;
mat[10] += mat[9] * mat[6] * di;
mat[11] += mat[9] * mat[7] * di;
mat[12] += mat[13] * mat[4] * di;
mat[14] += mat[13] * mat[6] * di;
mat[15] += mat[13] * mat[7] * di;
di = mat[10];
mat[10] = d = 1.0f / di;
mat[2] *= -d;
mat[6] *= -d;
mat[14] *= -d;
mat[8] *= d;
mat[9] *= d;
mat[11] *= d;
mat[0] += mat[2] * mat[8] * di;
mat[1] += mat[2] * mat[9] * di;
mat[3] += mat[2] * mat[11] * di;
mat[4] += mat[6] * mat[8] * di;
mat[5] += mat[6] * mat[9] * di;
mat[7] += mat[6] * mat[11] * di;
mat[12] += mat[14] * mat[8] * di;
mat[13] += mat[14] * mat[9] * di;
mat[15] += mat[14] * mat[11] * di;
di = mat[15];
mat[15] = d = 1.0f / di;
mat[3] *= -d;
mat[7] *= -d;
mat[11] *= -d;
mat[12] *= d;
mat[13] *= d;
mat[14] *= d;
mat[0] += mat[3] * mat[12] * di;
mat[1] += mat[3] * mat[13] * di;
mat[2] += mat[3] * mat[14] * di;
mat[4] += mat[7] * mat[12] * di;
mat[5] += mat[7] * mat[13] * di;
mat[6] += mat[7] * mat[14] * di;
mat[8] += mat[11] * mat[12] * di;
mat[9] += mat[11] * mat[13] * di;
mat[10] += mat[11] * mat[14] * di;
}

void PIII_inverse_4x4(float* src)
{
#if defined (MSCVER)
sf_m128 minor0, minor1, minor2, minor3;
sf_m128 row0, row1, row2, row3;
sf_m128 det, tmp1;
tmp1 = _mm_loadh_pi(_mm_loadl_pi(tmp1, (__m64*)(src) ), (__m64*)(src+ 4));
row1 = _mm_loadh_pi(_mm_loadl_pi(row1, (__m64*)(src+8)), (__m64*)(src+12));
row0 = _mm_shuffle_ps(tmp1, row1, 0x88);
row1 = _mm_shuffle_ps(row1, tmp1, 0xDD);
tmp1 = _mm_loadh_pi(_mm_loadl_pi(tmp1, (__m64*)(src+ 2)), (__m64*)(src+ 6));
row3 = _mm_loadh_pi(_mm_loadl_pi(row3, (__m64*)(src+10)), (__m64*)(src+14));
row2 = _mm_shuffle_ps(tmp1, row3, 0x88);
row3 = _mm_shuffle_ps(row3, tmp1, 0xDD);

// -----------------------------------------------

tmp1 = _mm_mul_ps(row2, row3);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
minor0 = _mm_mul_ps(row1, tmp1);
minor1 = _mm_mul_ps(row0, tmp1);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor0 = _mm_sub_ps(_mm_mul_ps(row1, tmp1), minor0);
minor1 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), minor1);
minor1 = _mm_shuffle_ps(minor1, minor1, 0x4E);

// -----------------------------------------------

tmp1 = _mm_mul_ps(row1, row2);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
minor0 = _mm_add_ps(_mm_mul_ps(row3, tmp1), minor0);
minor3 = _mm_mul_ps(row0, tmp1);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor0 = _mm_sub_ps(minor0, _mm_mul_ps(row3, tmp1));
minor3 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), minor3);
minor3 = _mm_shuffle_ps(minor3, minor3, 0x4E);

// -----------------------------------------------

tmp1 = _mm_mul_ps(_mm_shuffle_ps(row1, row1, 0x4E), row3);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
row2 = _mm_shuffle_ps(row2, row2, 0x4E);
minor0 = _mm_add_ps(_mm_mul_ps(row2, tmp1), minor0);
minor2 = _mm_mul_ps(row0, tmp1);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor0 = _mm_sub_ps(minor0, _mm_mul_ps(row2, tmp1));
minor2 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), minor2);
minor2 = _mm_shuffle_ps(minor2, minor2, 0x4E);

// -----------------------------------------------

tmp1 = _mm_mul_ps(row0, row1);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
minor2 = _mm_add_ps(_mm_mul_ps(row3, tmp1), minor2);
minor3 = _mm_sub_ps(_mm_mul_ps(row2, tmp1), minor3);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor2 = _mm_sub_ps(_mm_mul_ps(row3, tmp1), minor2);
minor3 = _mm_sub_ps(minor3, _mm_mul_ps(row2, tmp1));

// -----------------------------------------------

tmp1 = _mm_mul_ps(row0, row3);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
minor1 = _mm_sub_ps(minor1, _mm_mul_ps(row2, tmp1));
minor2 = _mm_add_ps(_mm_mul_ps(row1, tmp1), minor2);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor1 = _mm_add_ps(_mm_mul_ps(row2, tmp1), minor1);
minor2 = _mm_sub_ps(minor2, _mm_mul_ps(row1, tmp1));

// -----------------------------------------------

tmp1 = _mm_mul_ps(row0, row2);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
minor1 = _mm_add_ps(_mm_mul_ps(row3, tmp1), minor1);
minor3 = _mm_sub_ps(minor3, _mm_mul_ps(row1, tmp1));
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor1 = _mm_sub_ps(minor1, _mm_mul_ps(row3, tmp1));
minor3 = _mm_add_ps(_mm_mul_ps(row1, tmp1), minor3);

// -----------------------------------------------

det = _mm_mul_ps(row0, minor0);
det = _mm_add_ps(_mm_shuffle_ps(det, det, 0x4E), det);
det = _mm_add_ss(_mm_shuffle_ps(det, det, 0xB1), det);
tmp1 = _mm_rcp_ss(det);
det = _mm_sub_ss(_mm_add_ss(tmp1, tmp1), _mm_mul_ss(det, _mm_mul_ss(tmp1, tmp1)));
det = _mm_shuffle_ps(det, det, 0x00);
minor0 = _mm_mul_ps(det, minor0);
_mm_storel_pi((__m64*)(src), minor0);
_mm_storeh_pi((__m64*)(src+2), minor0);
minor1 = _mm_mul_ps(det, minor1);
_mm_storel_pi((__m64*)(src+4), minor1);
_mm_storeh_pi((__m64*)(src+6), minor1);
minor2 = _mm_mul_ps(det, minor2);
_mm_storel_pi((__m64*)(src+ 8), minor2);
_mm_storeh_pi((__m64*)(src+10), minor2);
minor3 = _mm_mul_ps(det, minor3);
_mm_storel_pi((__m64*)(src+12), minor3);
_mm_storeh_pi((__m64*)(src+14), minor3);
#endif
}

void PII_inverse_6x6(float* mat)
{
float d, di;
di = mat[0];
mat[0] = d = 1.0f / di;
mat[6] *= -d;
mat[12] *= -d;
mat[18] *= -d;
mat[24] *= -d;
mat[30] *= -d;
mat[1] *= d;
mat[2] *= d;
mat[3] *= d;
mat[4] *= d;
mat[5] *= d;
mat[7] += mat[6] * mat[1] * di;
mat[8] += mat[6] * mat[2] * di;
mat[9] += mat[6] * mat[3] * di;
mat[10] += mat[6] * mat[4] * di;
mat[11] += mat[6] * mat[5] * di;
mat[13] += mat[12] * mat[1] * di;
mat[14] += mat[12] * mat[2] * di;
mat[15] += mat[12] * mat[3] * di;
mat[16] += mat[12] * mat[4] * di;
mat[17] += mat[12] * mat[5] * di;
mat[19] += mat[18] * mat[1] * di;
mat[20] += mat[18] * mat[2] * di;
mat[21] += mat[18] * mat[3] * di;
mat[22] += mat[18] * mat[4] * di;
mat[23] += mat[18] * mat[5] * di;
mat[25] += mat[24] * mat[1] * di;
mat[26] += mat[24] * mat[2] * di;
mat[27] += mat[24] * mat[3] * di;
mat[28] += mat[24] * mat[4] * di;
mat[29] += mat[24] * mat[5] * di;
mat[31] += mat[30] * mat[1] * di;
mat[32] += mat[30] * mat[2] * di;
mat[33] += mat[30] * mat[3] * di;
mat[34] += mat[30] * mat[4] * di;
mat[35] += mat[30] * mat[5] * di;
di = mat[7];
mat[7] = d = 1.0f / di;
mat[1] *= -d;
mat[13] *= -d;
mat[19] *= -d;
mat[25] *= -d;
mat[31] *= -d;
mat[6] *= d;
mat[8] *= d;
mat[9] *= d;
mat[10] *= d;
mat[11] *= d;
mat[0] += mat[1] * mat[6] * di;
mat[2] += mat[1] * mat[8] * di;
mat[3] += mat[1] * mat[9] * di;
mat[4] += mat[1] * mat[10] * di;
mat[5] += mat[1] * mat[11] * di;
mat[12] += mat[13] * mat[6] * di;
mat[14] += mat[13] * mat[8] * di;
mat[15] += mat[13] * mat[9] * di;
mat[16] += mat[13] * mat[10] * di;
mat[17] += mat[13] * mat[11] * di;
mat[18] += mat[19] * mat[6] * di;
mat[20] += mat[19] * mat[8] * di;
mat[21] += mat[19] * mat[9] * di;
mat[22] += mat[19] * mat[10] * di;
mat[23] += mat[19] * mat[11] * di;
mat[24] += mat[25] * mat[6] * di;
mat[26] += mat[25] * mat[8] * di;
mat[27] += mat[25] * mat[9] * di;
mat[28] += mat[25] * mat[10] * di;
mat[29] += mat[25] * mat[11] * di;
mat[30] += mat[31] * mat[6] * di;
mat[32] += mat[31] * mat[8] * di;
mat[33] += mat[31] * mat[9] * di;

mat[34] += mat[31] * mat[10] * di;
mat[35] += mat[31] * mat[11] * di;
di = mat[14];
mat[14] = d = 1.0f / di;
mat[2] *= -d;
mat[8] *= -d;
mat[20] *= -d;
mat[26] *= -d;
mat[32] *= -d;
mat[12] *= d;
mat[13] *= d;
mat[15] *= d;
mat[16] *= d;
mat[17] *= d;
mat[0] += mat[2] * mat[12] * di;
mat[1] += mat[2] * mat[13] * di;
mat[3] += mat[2] * mat[15] * di;
mat[4] += mat[2] * mat[16] * di;
mat[5] += mat[2] * mat[17] * di;
mat[6] += mat[8] * mat[12] * di;
mat[7] += mat[8] * mat[13] * di;
mat[9] += mat[8] * mat[15] * di;
mat[10] += mat[8] * mat[16] * di;
mat[11] += mat[8] * mat[17] * di;
mat[18] += mat[20] * mat[12] * di;
mat[19] += mat[20] * mat[13] * di;
mat[21] += mat[20] * mat[15] * di;
mat[22] += mat[20] * mat[16] * di;
mat[23] += mat[20] * mat[17] * di;
mat[24] += mat[26] * mat[12] * di;
mat[25] += mat[26] * mat[13] * di;
mat[27] += mat[26] * mat[15] * di;
mat[28] += mat[26] * mat[16] * di;
mat[29] += mat[26] * mat[17] * di;
mat[30] += mat[32] * mat[12] * di;
mat[31] += mat[32] * mat[13] * di;
mat[33] += mat[32] * mat[15] * di;
mat[34] += mat[32] * mat[16] * di;
mat[35] += mat[32] * mat[17] * di;
di = mat[21];
mat[21] = d = 1.0f / di;
mat[3] *= -d;
mat[9] *= -d;
mat[15] *= -d;
mat[27] *= -d;
mat[33] *= -d;
mat[18] *= d;
mat[19] *= d;
mat[20] *= d;
mat[22] *= d;
mat[23] *= d;
mat[0] += mat[3] * mat[18] * di;
mat[1] += mat[3] * mat[19] * di;
mat[2] += mat[3] * mat[20] * di;
mat[4] += mat[3] * mat[22] * di;
mat[5] += mat[3] * mat[23] * di;
mat[6] += mat[9] * mat[18] * di;
mat[7] += mat[9] * mat[19] * di;
mat[8] += mat[9] * mat[20] * di;
mat[10] += mat[9] * mat[22] * di;
mat[11] += mat[9] * mat[23] * di;
mat[12] += mat[15] * mat[18] * di;
mat[13] += mat[15] * mat[19] * di;
mat[14] += mat[15] * mat[20] * di;
mat[16] += mat[15] * mat[22] * di;
mat[17] += mat[15] * mat[23] * di;
mat[24] += mat[27] * mat[18] * di;
mat[25] += mat[27] * mat[19] * di;
mat[26] += mat[27] * mat[20] * di;
mat[28] += mat[27] * mat[22] * di;
mat[29] += mat[27] * mat[23] * di;
mat[30] += mat[33] * mat[18] * di;
mat[31] += mat[33] * mat[19] * di;
mat[32] += mat[33] * mat[20] * di;
mat[34] += mat[33] * mat[22] * di;
mat[35] += mat[33] * mat[23] * di;
di = mat[28];
mat[28] = d = 1.0f / di;
mat[4] *= -d;
mat[10] *= -d;
mat[16] *= -d;
mat[22] *= -d;
mat[34] *= -d;
mat[24] *= d;
mat[25] *= d;
mat[26] *= d;
mat[27] *= d;
mat[29] *= d;
mat[0] += mat[4] * mat[24] * di;
mat[1] += mat[4] * mat[25] * di;
mat[2] += mat[4] * mat[26] * di;
mat[3] += mat[4] * mat[27] * di;
mat[5] += mat[4] * mat[29] * di;
mat[6] += mat[10] * mat[24] * di;
mat[7] += mat[10] * mat[25] * di;
mat[8] += mat[10] * mat[26] * di;
mat[9] += mat[10] * mat[27] * di;
mat[11] += mat[10] * mat[29] * di;
mat[12] += mat[16] * mat[24] * di;
mat[13] += mat[16] * mat[25] * di;
mat[14] += mat[16] * mat[26] * di;
mat[15] += mat[16] * mat[27] * di;
mat[17] += mat[16] * mat[29] * di;
mat[18] += mat[22] * mat[24] * di;
mat[19] += mat[22] * mat[25] * di;
mat[20] += mat[22] * mat[26] * di;
mat[21] += mat[22] * mat[27] * di;
mat[23] += mat[22] * mat[29] * di;
mat[30] += mat[34] * mat[24] * di;
mat[31] += mat[34] * mat[25] * di;
mat[32] += mat[34] * mat[26] * di;
mat[33] += mat[34] * mat[27] * di;
mat[35] += mat[34] * mat[29] * di;
di = mat[35];
mat[35] = d = 1.0f / di;
mat[5] *= -d;
mat[11] *= -d;
mat[17] *= -d;
mat[23] *= -d;
mat[29] *= -d;
mat[30] *= d;
mat[31] *= d;
mat[32] *= d;
mat[33] *= d;
mat[34] *= d;
mat[0] += mat[5] * mat[30] * di;
mat[1] += mat[5] * mat[31] * di;
mat[2] += mat[5] * mat[32] * di;
mat[3] += mat[5] * mat[33] * di;
mat[4] += mat[5] * mat[34] * di;
mat[6] += mat[11] * mat[30] * di;
mat[7] += mat[11] * mat[31] * di;
mat[8] += mat[11] * mat[32] * di;
mat[9] += mat[11] * mat[33] * di;
mat[10] += mat[11] * mat[34] * di;
mat[12] += mat[17] * mat[30] * di;
mat[13] += mat[17] * mat[31] * di;
mat[14] += mat[17] * mat[32] * di;
mat[15] += mat[17] * mat[33] * di;
mat[16] += mat[17] * mat[34] * di;
mat[18] += mat[23] * mat[30] * di;
mat[19] += mat[23] * mat[31] * di;
mat[20] += mat[23] * mat[32] * di;
mat[21] += mat[23] * mat[33] * di;
mat[22] += mat[23] * mat[34] * di;
mat[24] += mat[29] * mat[30] * di;
mat[25] += mat[29] * mat[31] * di;
mat[26] += mat[29] * mat[32] * di;
mat[27] += mat[29] * mat[33] * di;
mat[28] += mat[29] * mat[34] * di;
}

void PIII_InverseG_6x6(float *src)
{
#if defined(MSCVER)
#define EPSILON 1e-8
#define REAL_ZERO(x) (fabs(x) < EPSILON ? 1:0)
sf_m128 minor0, minor1, minor2, minor3;
sf_m128 det, tmp1, tmp2, tmp3, mask, index;
sf_m128 b[6];
sf_m128 row[6];
static const unsigned long minus_hex = 0x80000000;
static const sf_m128 minus = _mm_set_ps1(*(float*)&minus_hex);
static const sf_m128 e = _mm_set_ps(1.0f, 0.0f, 0.0f, 1.0f);
static const sf_m128 epsilon = _mm_set_ss(EPSILON);
float max, f;
int i, j, n1, n2, k, mask1, mask2, mask3;
// Loading matrixes: 4x2 to row[0], row[1] and 4x4 to row[2]...row[5].

tmp1 = _mm_loadh_pi(_mm_loadl_pi(tmp1, (__m64*)(&src[12])), (__m64*)(&src[18]));
tmp2 = _mm_loadh_pi(_mm_loadl_pi(tmp2, (__m64*)(&src[24])), (__m64*)(&src[30]));
row[0] = _mm_shuffle_ps(tmp1, tmp2, 0x88);
row[1] = _mm_shuffle_ps(tmp1, tmp2, 0xDD);
tmp1 = _mm_loadh_pi(_mm_loadl_pi(tmp1, (__m64*)(&src[14])), (__m64*)(&src[20]));
tmp2 = _mm_loadh_pi(_mm_loadl_pi(tmp2, (__m64*)(&src[26])), (__m64*)(&src[32]));
row[2] = _mm_shuffle_ps(tmp1, tmp2, 0x88);
row[3] = _mm_shuffle_ps(tmp1, tmp2, 0xDD);
tmp1 = _mm_loadh_pi(_mm_loadl_pi(tmp1, (__m64*)(&src[16])), (__m64*)(&src[22]));
tmp2 = _mm_loadh_pi(_mm_loadl_pi(tmp2, (__m64*)(&src[28])), (__m64*)(&src[34]));
row[4] = _mm_shuffle_ps(tmp1, tmp2, 0x88);
row[5] = _mm_shuffle_ps(tmp1, tmp2, 0xDD);

// Finding the max(|src[0]|, |src[1]|, ..., |src[5]|).

tmp1 = _mm_loadh_pi(_mm_load_ss(&src[2]), (__m64*)&src[0]);
tmp2 = _mm_loadh_pi(_mm_load_ss(&src[3]), (__m64*)&src[4]);
tmp1 = _mm_andnot_ps(minus, tmp1);
tmp2 = _mm_andnot_ps(minus, tmp2);
tmp3 = _mm_max_ps(tmp1, tmp2);
tmp3 = _mm_max_ps(tmp3, _mm_shuffle_ps(tmp3, tmp3, _MM_SHUFFLE(3, 2, 3, 2)));
tmp3 = _mm_max_ss(tmp3, _mm_shuffle_ps(tmp3, tmp3, _MM_SHUFFLE(1, 1, 1, 1)));
tmp3 = _mm_shuffle_ps(tmp3, tmp3, _MM_SHUFFLE(0, 0, 0, 0));
mask1 = _mm_movemask_ps(_mm_cmpeq_ps(tmp1, tmp3));
mask1 |= _mm_movemask_ps(_mm_cmpeq_ps(tmp2, tmp3))<<4;
mask2 = mask1 & 0x98;
mask2 = mask2 - (mask2 << 1);
n1 = ((unsigned int)mask2) >> 31;
n1 |= ((mask1 & 0x11) != 0) << 1;
mask2 = mask1 & 0xC0;
mask2 = mask2 - (mask2 << 1);
n1 |= (((unsigned int)mask2) >> 29) & 4;
if(REAL_ZERO(src[n1]))
return;

// The first Gauss iteration.

tmp1 = row[n1];
row[n1] = row[0];
row[0] = tmp1;
tmp2 = _mm_load_ss(&src[n1]);
src[n1] = src[0];
f = src[n1+6];
src[n1+6] = src[6];
src[6] = f;
tmp1 = _mm_rcp_ss(tmp2);
tmp2 = _mm_sub_ss(_mm_add_ss(tmp1, tmp1), _mm_mul_ss(tmp2, _mm_mul_ss(tmp1, tmp1)));
_mm_store_ss(&src[0], tmp2);
tmp2 = _mm_shuffle_ps(tmp2, tmp2, 0x00);
row[0] = _mm_mul_ps(row[0], tmp2);
tmp1 = _mm_load_ss(&src[1]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
row[1] = _mm_sub_ps(row[1], _mm_mul_ps(row[0], tmp1));
tmp1 = _mm_load_ss(&src[2]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
row[2] = _mm_sub_ps(row[2], _mm_mul_ps(row[0], tmp1));
tmp1 = _mm_load_ss(&src[3]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
row[3] = _mm_sub_ps(row[3], _mm_mul_ps(row[0], tmp1));
tmp1 = _mm_load_ss(&src[4]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
row[4] = _mm_sub_ps(row[4], _mm_mul_ps(row[0], tmp1));
tmp1 = _mm_load_ss(&src[5]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
row[5] = _mm_sub_ps(row[5], _mm_mul_ps(row[0], tmp1));
tmp3 = _mm_load_ss(&src[6]);
tmp3 = _mm_mul_ss(tmp3, tmp2);
_mm_store_ss(&src[6], tmp3);
tmp1 = _mm_load_ss(&src[1]);
tmp2 = _mm_load_ss(&src[7]);
tmp2 = _mm_sub_ss(tmp2, _mm_mul_ss(tmp1, tmp3));
_mm_store_ss(&src[7], tmp2);
tmp3 = _mm_shuffle_ps(tmp3, tmp3, 0x00);
tmp1 = _mm_loadh_pi(_mm_loadl_pi(tmp1, (__m64*)&src[2]), (__m64*)&src[ 4]);
tmp2 = _mm_loadh_pi(_mm_loadl_pi(tmp2, (__m64*)&src[8]), (__m64*)&src[10]);
tmp2 = _mm_sub_ps(tmp2, _mm_mul_ps(tmp1, tmp3));
_mm_storel_pi((__m64*)&src[ 8], tmp2);
_mm_storeh_pi((__m64*)&src[10], tmp2);

// Finding the max(src[7], src[8], ..., src[11]).

tmp1 = _mm_loadh_pi(_mm_load_ss(&src[7]), (__m64*)&src[10]);
tmp2 = _mm_loadl_pi(tmp2, (__m64*)&src[8]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, _MM_SHUFFLE(0,3,2,2));
tmp1 = _mm_andnot_ps(minus, tmp1);
tmp2 = _mm_andnot_ps(minus, tmp2);
tmp3 = _mm_max_ps(tmp1, tmp2);
tmp3 = _mm_max_ps(tmp3, _mm_shuffle_ps(tmp1, tmp1, _MM_SHUFFLE(0,0,3,2)));
tmp3 = _mm_max_ss(tmp3, _mm_shuffle_ps(tmp3, tmp3, _MM_SHUFFLE(1,1,1,1)));
tmp3 = _mm_shuffle_ps(tmp3, tmp3, _MM_SHUFFLE(0,0,0,0));
mask1 = _mm_movemask_ps(_mm_cmpeq_ps(tmp2, tmp3));
mask2 = _mm_movemask_ps(_mm_cmpeq_ps(tmp1, tmp3));
n2 = ((mask1 & 3) | (mask2 & 7)) + 7;
if(REAL_ZERO(src[n2]))
return;

// The second Gauss iteration.

tmp2 = _mm_load_ss(&src[n2]);
src[n2] = src[7];
n2 -= 6;
tmp1 = row[n2];
row[n2] = row[1];
row[1] = tmp1;
f = src[n2];
src[n2] = src[1];
src[1] = f;
//if(n2==n1) n2 = 0;

n2 *= (n1!=n2);
tmp1 = _mm_rcp_ss(tmp2);
tmp2 = _mm_sub_ss(_mm_add_ss(tmp1, tmp1), _mm_mul_ss(tmp2, _mm_mul_ss(tmp1, tmp1)));
_mm_store_ss(&src[7], tmp2);
tmp2 = _mm_shuffle_ps(tmp2, tmp2, 0x00);
row[1] = _mm_mul_ps(row[1], tmp2);
tmp1 = _mm_load_ss(&src[6]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
row[0] = _mm_sub_ps(row[0], _mm_mul_ps(row[1], tmp1));
tmp1 = _mm_load_ss(&src[8]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
row[2] = _mm_sub_ps(row[2], _mm_mul_ps(row[1], tmp1));
tmp1 = _mm_load_ss(&src[9]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
row[3] = _mm_sub_ps(row[3], _mm_mul_ps(row[1], tmp1));
tmp1 = _mm_load_ss(&src[10]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
row[4] = _mm_sub_ps(row[4], _mm_mul_ps(row[1], tmp1));
tmp1 = _mm_load_ss(&src[11]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
row[5] = _mm_sub_ps(row[5], _mm_mul_ps(row[1], tmp1));
row[0] = _mm_xor_ps(row[0], minus);
row[1] = _mm_xor_ps(row[1], minus);

// Inverting the matrix 4x4 by the Kramers method.

row[3] = _mm_shuffle_ps(row[3], row[3], 0x4E);
row[5] = _mm_shuffle_ps(row[5], row[5], 0x4E);
tmp2 = _mm_mul_ps(row[4], row[5]);
tmp1 = _mm_shuffle_ps(tmp2, tmp2, 0xB1);
minor0 = _mm_mul_ps(row[3], tmp1);
minor1 = _mm_mul_ps(row[2], tmp1);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor0 = _mm_sub_ps(_mm_mul_ps(row[3], tmp1), minor0);
minor1 = _mm_sub_ps(_mm_mul_ps(row[2], tmp1), minor1);
minor1 = _mm_shuffle_ps(minor1, minor1, 0x4E);

// -----------------------------------------------

tmp2 = _mm_mul_ps(row[3], row[4]);
tmp1 = _mm_shuffle_ps(tmp2, tmp2, 0xB1);
minor0 = _mm_add_ps(_mm_mul_ps(row[5], tmp1), minor0);
minor3 = _mm_mul_ps(row[2], tmp1);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor0 = _mm_sub_ps(minor0, _mm_mul_ps(row[5], tmp1));
minor3 = _mm_sub_ps(_mm_mul_ps(row[2], tmp1), minor3);
minor3 = _mm_shuffle_ps(minor3, minor3, 0x4E);

// -----------------------------------------------

tmp2 = _mm_mul_ps(_mm_shuffle_ps(row[3], row[3], 0x4E), row[5]);
tmp1 = _mm_shuffle_ps(tmp2, tmp2, 0xB1);
row[4] = _mm_shuffle_ps(row[4], row[4], 0x4E);
minor0 = _mm_add_ps(_mm_mul_ps(row[4], tmp1), minor0);
minor2 = _mm_mul_ps(row[2], tmp1);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor0 = _mm_sub_ps(minor0, _mm_mul_ps(row[4], tmp1));
minor2 = _mm_sub_ps(_mm_mul_ps(row[2], tmp1), minor2);
minor2 = _mm_shuffle_ps(minor2, minor2, 0x4E);

// -----------------------------------------------

tmp2 = _mm_mul_ps(row[2], row[3]);
tmp1 = _mm_shuffle_ps(tmp2, tmp2, 0xB1);
minor2 = _mm_add_ps(_mm_mul_ps(row[5], tmp1), minor2);
minor3 = _mm_sub_ps(_mm_mul_ps(row[4], tmp1), minor3);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor2 = _mm_sub_ps(_mm_mul_ps(row[5], tmp1), minor2);
minor3 = _mm_sub_ps(minor3, _mm_mul_ps(row[4], tmp1));

// -----------------------------------------------

tmp2 = _mm_mul_ps(row[2], row[5]);
tmp1 = _mm_shuffle_ps(tmp2, tmp2, 0xB1);
minor1 = _mm_sub_ps(minor1, _mm_mul_ps(row[4], tmp1));
minor2 = _mm_add_ps(_mm_mul_ps(row[3], tmp1), minor2);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor1 = _mm_add_ps(_mm_mul_ps(row[4], tmp1), minor1);
minor2 = _mm_sub_ps(minor2, _mm_mul_ps(row[3], tmp1));

// -----------------------------------------------

tmp2 = _mm_mul_ps(row[2], row[4]);
tmp1 = _mm_shuffle_ps(tmp2, tmp2, 0xB1);
minor1 = _mm_add_ps(_mm_mul_ps(row[5], tmp1), minor1);
minor3 = _mm_sub_ps(minor3, _mm_mul_ps(row[3], tmp1));
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor1 = _mm_sub_ps(minor1, _mm_mul_ps(row[5], tmp1));
minor3 = _mm_add_ps(_mm_mul_ps(row[3], tmp1), minor3);

// -----------------------------------------------

det = _mm_mul_ps(row[2], minor0);
det = _mm_add_ps(_mm_shuffle_ps(det, det, 0x4E), det);
det = _mm_add_ss(_mm_shuffle_ps(det, det, 0xB1), det);
if(_mm_movemask_ps(_mm_cmplt_ss(_mm_andnot_ps(minus, det), epsilon)) & 1)
return;
tmp1 = _mm_rcp_ss(det);
det = _mm_sub_ss(_mm_add_ss(tmp1, tmp1), _mm_mul_ss(det, _mm_mul_ss(tmp1, tmp1)));
det = _mm_shuffle_ps(det, det, 0x00);
row[2] = _mm_mul_ps(det, minor0);
row[3] = _mm_mul_ps(det, minor1);

///////////////////////////////////////

b[0] = _mm_unpacklo_ps(row[0], row[1]);
b[2] = _mm_unpackhi_ps(row[0], row[1]);
row[4] = _mm_mul_ps(det, minor2);

b[1] = _mm_shuffle_ps(b[0], b[2], 0x4E);
row[5] = _mm_mul_ps(det, minor3);
b[3] = _mm_shuffle_ps(b[2], b[0], 0x4E);
tmp1 = _mm_shuffle_ps(row[2], row[3], 0x50);
tmp2 = _mm_mul_ps(b[0], tmp1);
tmp1 = _mm_shuffle_ps(row[2], row[3], 0xA5);
tmp2 = _mm_add_ps(tmp2, _mm_mul_ps(b[1], tmp1));
tmp1 = _mm_shuffle_ps(row[2], row[3], 0xFA);
tmp2 = _mm_add_ps(tmp2, _mm_mul_ps(b[2], tmp1));
tmp1 = _mm_shuffle_ps(row[2], row[3], 0x0F);
row[0] = _mm_add_ps(tmp2, _mm_mul_ps(b[3], tmp1));
tmp1 = _mm_shuffle_ps(row[4], row[5], 0x50);
tmp2 = _mm_mul_ps(b[0], tmp1);
tmp1 = _mm_shuffle_ps(row[4], row[5], 0xA5);
tmp2 = _mm_add_ps(tmp2, _mm_mul_ps(b[1], tmp1));
tmp1 = _mm_shuffle_ps(row[4], row[5], 0xFA);
tmp2 = _mm_add_ps(tmp2, _mm_mul_ps(b[2], tmp1));
tmp1 = _mm_shuffle_ps(row[4], row[5], 0x0F);
row[1] = _mm_add_ps(tmp2, _mm_mul_ps(b[3], tmp1));
b[2] = _mm_shuffle_ps(row[0], row[0], 0x44);
b[3] = _mm_shuffle_ps(row[0], row[0], 0xEE);
b[4] = _mm_shuffle_ps(row[1], row[1], 0x44);
b[5] = _mm_shuffle_ps(row[1], row[1], 0xEE);

// Calculating row number n2

tmp1 = _mm_load_ss(&src[8]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
b [1] = _mm_sub_ps(_mm_shuffle_ps(e, e, 0x4E), _mm_mul_ps(b[2], tmp1));
row[1] = _mm_xor_ps(_mm_mul_ps(row[2], tmp1), minus);
tmp1 = _mm_load_ss(&src[9]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
b [1] = _mm_sub_ps(b [1], _mm_mul_ps(b [3], tmp1));
row[1] = _mm_sub_ps(row[1], _mm_mul_ps(row[3], tmp1));
tmp1 = _mm_load_ss(&src[10]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
b [1] = _mm_sub_ps(b [1], _mm_mul_ps(b [4], tmp1));
row[1] = _mm_sub_ps(row[1], _mm_mul_ps(row[4], tmp1));
tmp1 = _mm_load_ss(&src[11]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
b [1] = _mm_sub_ps(b [1], _mm_mul_ps(b [5], tmp1));
row[1] = _mm_sub_ps(row[1], _mm_mul_ps(row[5], tmp1));
tmp1 = _mm_load_ss(&src[6]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
b [1] = _mm_sub_ps(b[1], _mm_mul_ps(e, tmp1));
tmp2 = _mm_load_ss(&src[7]);
tmp2 = _mm_shuffle_ps(tmp2, tmp2, 0x00);
b [1] = _mm_mul_ps(b [1], tmp2);
row[1] = _mm_mul_ps(row[1], tmp2);

// Calculating row number n1

tmp1 = _mm_load_ss(&src[1]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
b [0] = _mm_sub_ps(e, _mm_mul_ps(b[1], tmp1));
row[0] = _mm_xor_ps(_mm_mul_ps(row[1], tmp1), minus);
tmp1 = _mm_load_ss(&src[2]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
b [0] = _mm_sub_ps(b [0], _mm_mul_ps(b [2], tmp1));
row[0] = _mm_sub_ps(row[0], _mm_mul_ps(row[2], tmp1));
tmp1 = _mm_load_ss(&src[3]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
b [0] = _mm_sub_ps(b [0], _mm_mul_ps(b [3], tmp1));
row[0] = _mm_sub_ps(row[0], _mm_mul_ps(row[3], tmp1));
tmp1 = _mm_load_ss(&src[4]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
b [0] = _mm_sub_ps(b [0], _mm_mul_ps(b [4], tmp1));
row[0] = _mm_sub_ps(row[0], _mm_mul_ps(row[4], tmp1));
tmp1 = _mm_load_ss(&src[5]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
b [0] = _mm_sub_ps(b [0], _mm_mul_ps(b [5], tmp1));
row[0] = _mm_sub_ps(row[0], _mm_mul_ps(row[5], tmp1));
tmp2 = _mm_load_ss(&src[0]);
tmp2 = _mm_shuffle_ps(tmp2, tmp2, 0x00);
b [0] = _mm_mul_ps(b [0], tmp2);
row[0] = _mm_mul_ps(row[0], tmp2);
n2 = (n2==0)*(n1-n2)+n2;
tmp1 = row[ 1]; row[ 1] = row[n2]; row[n2] = tmp1;
tmp2 = b [ 1]; b [ 1] = b [n2]; b [n2] = tmp2;
tmp1 = row[ 0]; row[ 0] = row[n1]; row[n1] = tmp1;
tmp2 = b [ 0]; b [ 0] = b [n1]; b [n1] = tmp2;
_mm_storel_pi((__m64*)&src[ 0], b [0]);
_mm_storel_pi((__m64*)&src[ 2], row[0]);
_mm_storeh_pi((__m64*)&src[ 4], row[0]);
_mm_storel_pi((__m64*)&src[ 6], b [1]);
_mm_storel_pi((__m64*)&src[ 8], row[1]);
_mm_storeh_pi((__m64*)&src[10], row[1]);
_mm_storel_pi((__m64*)&src[12], b [2]);
_mm_storel_pi((__m64*)&src[14], row[2]);
_mm_storeh_pi((__m64*)&src[16], row[2]);
_mm_storel_pi((__m64*)&src[18], b [3]);
_mm_storel_pi((__m64*)&src[20], row[3]);
_mm_storeh_pi((__m64*)&src[22], row[3]);
_mm_storel_pi((__m64*)&src[24], b [4]);
_mm_storel_pi((__m64*)&src[26], row[4]);
_mm_storeh_pi((__m64*)&src[28], row[4]);
_mm_storel_pi((__m64*)&src[30], b [5]);
_mm_storel_pi((__m64*)&src[32], row[5]);
_mm_storeh_pi((__m64*)&src[34], row[5]);
#undef EPSILON
#undef REAL_ZERO
#endif

}

// PIII_InverseG_6x6

void PIII_InverseS_6x6(float *src)
{
#if defined (MSCVER)
#define EPSILON 1e-8
#define REAL_ZERO(x) (fabs(x) < EPSILON ? 1:0)

sf_m128 minor0, minor1, minor2, minor3;
sf_m128 det, tmp1, tmp2;
sf_m128 b0, b1, b2, b3;
sf_m128 row[6];
static const unsigned long minus_hex = 0x80000000;
static const sf_m128 minus = _mm_set_ps1(*(float*)&minus_hex);
static const sf_m128 zero = _mm_setzero_ps();
static const sf_m128 e = _mm_set_ps(1.0f,0.0f,0.0f,1.0f);
static const sf_m128 epsilon = _mm_set_ss(EPSILON);
static const sf_m128 epsilon1 = _mm_set_ss(-EPSILON);

// Loading matrixes: 4x2 to row[0], row[1] and 4x4 to row[2]...row[5].

tmp1 = _mm_loadh_pi(_mm_loadl_pi(tmp1, (__m64*)(&src[12])), (__m64*)(&src[18]));
tmp2 = _mm_loadh_pi(_mm_loadl_pi(tmp2, (__m64*)(&src[24])), (__m64*)(&src[30]));
row[0] = _mm_shuffle_ps(tmp1, tmp2, 0x88);
row[1] = _mm_shuffle_ps(tmp1, tmp2, 0xDD);
tmp1 = _mm_loadh_pi(_mm_loadl_pi(tmp1, (__m64*)(&src[14])), (__m64*)(&src[20]));
tmp2 = _mm_loadh_pi(_mm_loadl_pi(tmp2, (__m64*)(&src[26])), (__m64*)(&src[32]));
row[2] = _mm_shuffle_ps(tmp1, tmp2, 0x88);
row[3] = _mm_shuffle_ps(tmp1, tmp2, 0xDD);
tmp1 = _mm_loadh_pi(_mm_loadl_pi(tmp1, (__m64*)(&src[16])), (__m64*)(&src[22]));
tmp2 = _mm_loadh_pi(_mm_loadl_pi(tmp2, (__m64*)(&src[28])), (__m64*)(&src[34]));
row[4] = _mm_shuffle_ps(tmp1, tmp2, 0x88);
row[5] = _mm_shuffle_ps(tmp1, tmp2, 0xDD);

// ----------------

tmp2 = _mm_load_ss(&src[0]);
tmp1 = _mm_rcp_ss(tmp2);
tmp2 = _mm_sub_ss(_mm_add_ss(tmp1, tmp1), _mm_mul_ss(tmp2, _mm_mul_ss(tmp1, tmp1)));
_mm_store_ss(&src[0], tmp2);
tmp2 = _mm_shuffle_ps(tmp2, tmp2, 0x00);
row[0] = _mm_mul_ps(row[0], tmp2);
tmp1 = _mm_load_ss(&src[1]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
row[1] = _mm_sub_ps(row[1], _mm_mul_ps(row[0], tmp1));
tmp1 = _mm_load_ss(&src[2]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
row[2] = _mm_sub_ps(row[2], _mm_mul_ps(row[0], tmp1));
tmp1 = _mm_load_ss(&src[3]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
row[3] = _mm_sub_ps(row[3], _mm_mul_ps(row[0], tmp1));
tmp1 = _mm_load_ss(&src[4]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
row[4] = _mm_sub_ps(row[4], _mm_mul_ps(row[0], tmp1));
tmp1 = _mm_load_ss(&src[5]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
row[5] = _mm_sub_ps(row[5], _mm_mul_ps(row[0], tmp1));
b0 = _mm_load_ss(&src[6]);
b0 = _mm_mul_ss(b0, tmp2);
_mm_store_ss(&src[6], b0);
tmp1 = _mm_load_ss(&src[1]);
tmp2 = _mm_load_ss(&src[7]);
tmp2 = _mm_sub_ss(tmp2, _mm_mul_ss(tmp1, b0));
_mm_store_ss(&src[7], tmp2);
b0 = _mm_shuffle_ps(b0, b0, 0x00);
tmp1 = _mm_loadh_pi(_mm_loadl_pi(tmp1, (__m64*)&src[2]), (__m64*)&src[ 4]);
tmp2 = _mm_loadh_pi(_mm_loadl_pi(tmp2, (__m64*)&src[8]), (__m64*)&src[10]);
tmp2 = _mm_sub_ps(tmp2, _mm_mul_ps(tmp1, b0));
_mm_storel_pi((__m64*)&src[ 8], tmp2);
_mm_storeh_pi((__m64*)&src[10], tmp2);

// ----------------

tmp2 = _mm_load_ss(&src[7]);
tmp1 = _mm_rcp_ss(tmp2);
tmp2 = _mm_sub_ss(_mm_add_ss(tmp1, tmp1), _mm_mul_ss(tmp2, _mm_mul_ss(tmp1, tmp1)));
_mm_store_ss(&src[7], tmp2);
tmp2 = _mm_shuffle_ps(tmp2, tmp2, 0x00);
row[1] = _mm_mul_ps(row[1], tmp2);
row[0] = _mm_sub_ps(row[0], _mm_mul_ps(row[1], b0));
tmp1 = _mm_load_ss(&src[8]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
row[2] = _mm_sub_ps(row[2], _mm_mul_ps(row[1], tmp1));
tmp1 = _mm_load_ss(&src[9]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
row[3] = _mm_sub_ps(row[3], _mm_mul_ps(row[1], tmp1));
tmp1 = _mm_load_ss(&src[10]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
row[4] = _mm_sub_ps(row[4], _mm_mul_ps(row[1], tmp1));
tmp1 = _mm_load_ss(&src[11]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
row[5] = _mm_sub_ps(row[5], _mm_mul_ps(row[1], tmp1));
row[0] = _mm_xor_ps(row[0], minus);
row[1] = _mm_xor_ps(row[1], minus);
row[3] = _mm_shuffle_ps(row[3], row[3], 0x4E);
row[5] = _mm_shuffle_ps(row[5], row[5], 0x4E);

// Inverting the matrix 4x4 by the Kramers method.

tmp2 = _mm_mul_ps(row[4], row[5]);
tmp1 = _mm_shuffle_ps(tmp2, tmp2, 0xB1);
minor0 = _mm_mul_ps(row[3], tmp1);
minor1 = _mm_mul_ps(row[2], tmp1);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor0 = _mm_sub_ps(_mm_mul_ps(row[3], tmp1), minor0);
minor1 = _mm_sub_ps(_mm_mul_ps(row[2], tmp1), minor1);
minor1 = _mm_shuffle_ps(minor1, minor1, 0x4E);

// -----------------------------------------------

tmp2 = _mm_mul_ps(row[3], row[4]);
tmp1 = _mm_shuffle_ps(tmp2, tmp2, 0xB1);
minor0 = _mm_add_ps(_mm_mul_ps(row[5], tmp1), minor0);
minor3 = _mm_mul_ps(row[2], tmp1);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor0 = _mm_sub_ps(minor0, _mm_mul_ps(row[5], tmp1));
minor3 = _mm_sub_ps(_mm_mul_ps(row[2], tmp1), minor3);
minor3 = _mm_shuffle_ps(minor3, minor3, 0x4E);

// -----------------------------------------------

tmp2 = _mm_mul_ps(_mm_shuffle_ps(row[3], row[3], 0x4E), row[5]);
tmp1 = _mm_shuffle_ps(tmp2, tmp2, 0xB1);
row[4] = _mm_shuffle_ps(row[4], row[4], 0x4E);
minor0 = _mm_add_ps(_mm_mul_ps(row[4], tmp1), minor0);
minor2 = _mm_mul_ps(row[2], tmp1);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor0 = _mm_sub_ps(minor0, _mm_mul_ps(row[4], tmp1));
minor2 = _mm_sub_ps(_mm_mul_ps(row[2], tmp1), minor2);
minor2 = _mm_shuffle_ps(minor2, minor2, 0x4E);

// -----------------------------------------------

tmp2 = _mm_mul_ps(row[2], row[3]);
tmp1 = _mm_shuffle_ps(tmp2, tmp2, 0xB1);
minor2 = _mm_add_ps(_mm_mul_ps(row[5], tmp1), minor2);
minor3 = _mm_sub_ps(_mm_mul_ps(row[4], tmp1), minor3);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor2 = _mm_sub_ps(_mm_mul_ps(row[5], tmp1), minor2);
minor3 = _mm_sub_ps(minor3, _mm_mul_ps(row[4], tmp1));

// -----------------------------------------------

tmp2 = _mm_mul_ps(row[2], row[5]);
tmp1 = _mm_shuffle_ps(tmp2, tmp2, 0xB1);
minor1 = _mm_sub_ps(minor1, _mm_mul_ps(row[4], tmp1));
minor2 = _mm_add_ps(_mm_mul_ps(row[3], tmp1), minor2);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor1 = _mm_add_ps(_mm_mul_ps(row[4], tmp1), minor1);
minor2 = _mm_sub_ps(minor2, _mm_mul_ps(row[3], tmp1));

// -----------------------------------------------

tmp2 = _mm_mul_ps(row[2], row[4]);
tmp1 = _mm_shuffle_ps(tmp2, tmp2, 0xB1);
minor1 = _mm_add_ps(_mm_mul_ps(row[5], tmp1), minor1);
minor3 = _mm_sub_ps(minor3, _mm_mul_ps(row[3], tmp1));
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor1 = _mm_sub_ps(minor1, _mm_mul_ps(row[5], tmp1));
minor3 = _mm_add_ps(_mm_mul_ps(row[3], tmp1), minor3);

// -----------------------------------------------

det = _mm_mul_ps(row[2], minor0);
det = _mm_add_ps(_mm_shuffle_ps(det, det, 0x4E), det);
det = _mm_add_ss(_mm_shuffle_ps(det, det, 0xB1), det);
if(_mm_movemask_ps(_mm_and_ps(_mm_cmplt_ss(det, epsilon), _mm_cmpgt_ss(det, epsilon1))) & 1)
return;

tmp1 = _mm_rcp_ss(det);
det = _mm_sub_ss(_mm_add_ss(tmp1, tmp1), _mm_mul_ss(det, _mm_mul_ss(tmp1, tmp1)));
det = _mm_shuffle_ps(det, det, 0x00);
row[2] = _mm_mul_ps(det, minor0);
row[3] = _mm_mul_ps(det, minor1);
row[4] = _mm_mul_ps(det, minor2);
row[5] = _mm_mul_ps(det, minor3);
b0 = _mm_unpacklo_ps(row[0], row[1]);
b2 = _mm_unpackhi_ps(row[0], row[1]);
b1 = _mm_shuffle_ps(b0, b2, 0x4E);
b3 = _mm_shuffle_ps(b2, b0, 0x4E);
tmp1 = _mm_shuffle_ps(row[2], row[3], 0x50);
tmp2 = _mm_mul_ps(b0, tmp1);
tmp1 = _mm_shuffle_ps(row[2], row[3], 0xA5);
tmp2 = _mm_add_ps(tmp2, _mm_mul_ps(b1, tmp1));
tmp1 = _mm_shuffle_ps(row[2], row[3], 0xFA);
tmp2 = _mm_add_ps(tmp2, _mm_mul_ps(b2, tmp1));
tmp1 = _mm_shuffle_ps(row[2], row[3], 0x0F);
row[0] = _mm_add_ps(tmp2, _mm_mul_ps(b3, tmp1));
tmp1 = _mm_shuffle_ps(row[4], row[5], 0x50);
tmp2 = _mm_mul_ps(b0, tmp1);
tmp1 = _mm_shuffle_ps(row[4], row[5], 0xA5);
tmp2 = _mm_add_ps(tmp2, _mm_mul_ps(b1, tmp1));
tmp1 = _mm_shuffle_ps(row[4], row[5], 0xFA);
tmp2 = _mm_add_ps(tmp2, _mm_mul_ps(b2, tmp1));
tmp1 = _mm_shuffle_ps(row[4], row[5], 0x0F);
row[1] = _mm_add_ps(tmp2, _mm_mul_ps(b3, tmp1));

// Calculating row number 1

b0 = e;
tmp1 = _mm_load_ss(&src[8]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
b0 = _mm_sub_ps(b0, _mm_mul_ps(_mm_shuffle_ps(row[0], row[0], 0x4E), tmp1));
b1 = _mm_xor_ps(_mm_mul_ps(row[2], tmp1), minus);
tmp1 = _mm_load_ss(&src[9]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
b0 = _mm_sub_ps(b0, _mm_mul_ps(row[0], tmp1));
b1 = _mm_sub_ps(b1, _mm_mul_ps(row[3], tmp1));
tmp1 = _mm_load_ss(&src[10]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
b0 = _mm_sub_ps(b0, _mm_mul_ps(_mm_shuffle_ps(row[1], row[1], 0x4E), tmp1));
b1 = _mm_sub_ps(b1, _mm_mul_ps(row[4], tmp1));
tmp1 = _mm_load_ss(&src[11]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
b0 = _mm_sub_ps(b0, _mm_mul_ps(row[1], tmp1));
b1 = _mm_sub_ps(b1, _mm_mul_ps(row[5], tmp1));
tmp1 = _mm_load_ss(&src[6]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
b0 = _mm_sub_ps(b0, _mm_mul_ps(_mm_shuffle_ps(e, e, 0x4E), tmp1));
tmp2 = _mm_load_ss(&src[7]);
tmp2 = _mm_shuffle_ps(tmp2, tmp2, 0x00);
b0 = _mm_mul_ps(b0, tmp2);
b1 = _mm_mul_ps(b1, tmp2);
_mm_storeh_pi((__m64*)&src[ 6], b0);
_mm_storel_pi((__m64*)&src[ 8], b1);
_mm_storeh_pi((__m64*)&src[10], b1);

// Calculating row number 0

tmp1 = _mm_load_ss(&src[1]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
b0 = _mm_sub_ps(e, _mm_mul_ps(_mm_shuffle_ps(b0, b0, 0x4E), tmp1));
b1 = _mm_xor_ps(_mm_mul_ps(b1, tmp1), minus);
tmp1 = _mm_load_ss(&src[2]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
b0 = _mm_sub_ps(b0, _mm_mul_ps(row[0], tmp1));
b1 = _mm_sub_ps(b1, _mm_mul_ps(row[2], tmp1));
tmp1 = _mm_load_ss(&src[3]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
b0 = _mm_sub_ps(b0, _mm_mul_ps(_mm_shuffle_ps(row[0], row[0], 0x4E), tmp1));
b1 = _mm_sub_ps(b1, _mm_mul_ps(row[3], tmp1));
tmp1 = _mm_load_ss(&src[4]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
b0 = _mm_sub_ps(b0, _mm_mul_ps(row[1], tmp1));
b1 = _mm_sub_ps(b1, _mm_mul_ps(row[4], tmp1));
tmp1 = _mm_load_ss(&src[5]);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x00);
b0 = _mm_sub_ps(b0, _mm_mul_ps(_mm_shuffle_ps(row[1], row[1], 0x4E), tmp1));
b1 = _mm_sub_ps(b1, _mm_mul_ps(row[5], tmp1));
tmp2 = _mm_load_ss(&src[0]);
tmp2 = _mm_shuffle_ps(tmp2, tmp2, 0x00);
b0 = _mm_mul_ps(b0, tmp2);
b1 = _mm_mul_ps(b1, tmp2);
_mm_storel_pi((__m64*)&src[0], b0);
_mm_storel_pi((__m64*)&src[2], b1);
_mm_storeh_pi((__m64*)&src[4], b1);
_mm_storel_pi((__m64*)&src[12], row[0]);
_mm_storel_pi((__m64*)&src[14], row[2]);
_mm_storeh_pi((__m64*)&src[16], row[2]);
_mm_storeh_pi((__m64*)&src[18], row[0]);
_mm_storel_pi((__m64*)&src[20], row[3]);
_mm_storeh_pi((__m64*)&src[22], row[3]);
_mm_storel_pi((__m64*)&src[24], row[1]);
_mm_storel_pi((__m64*)&src[26], row[4]);
_mm_storeh_pi((__m64*)&src[28], row[4]);
_mm_storeh_pi((__m64*)&src[30], row[1]);
_mm_storel_pi((__m64*)&src[32], row[5]);
_mm_storeh_pi((__m64*)&src[34], row[5]);
#endif
} // PIII_InverseS_6x6



//===============================================================
//
//	CMat2D
//
//===============================================================

CMat2D mat2_zero( CVec2D( 0, 0 ), CVec2D( 0, 0 ) );
CMat2D mat2_identity( CVec2D( 1, 0 ), CVec2D( 0, 1 ) );

/*
============
CMat2D::InverseSelf
============
*/
bool CMat2D::inverseSelf() {
	// 2+4 = 6 multiplications
	//		 1 division
	double det, invDet, a;

	det = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];

	if ( CMath::fabs( det ) < MATRIX_INVERSE_EPSILON ) {
		return false;
	}

	invDet = 1.0f / det;

	a = mat[0][0];
	mat[0][0] =   mat[1][1] * invDet;
	mat[0][1] = - mat[0][1] * invDet;
	mat[1][0] = - mat[1][0] * invDet;
	mat[1][1] =   a * invDet;

	return true;
}

/*
============
CMat2D::InverseFastSelf
============
*/
bool CMat2D::inverseFastSelf() {
#if 1
	// 2+4 = 6 multiplications
	//		 1 division
	double det, invDet, a;

	det = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];

	if ( CMath::fabs( det ) < MATRIX_INVERSE_EPSILON ) {
		return false;
	}

	invDet = 1.0f / det;

	a = mat[0][0];
	mat[0][0] =   mat[1][1] * invDet;
	mat[0][1] = - mat[0][1] * invDet;
	mat[1][0] = - mat[1][0] * invDet;
	mat[1][1] =   a * invDet;

	return true;
#else
	// 2*4 = 8 multiplications
	//		 2 division
	float *mat = reinterpret_cast<float *>(this);
	double d, di;
	float s;

	di = mat[0];
	s = di;
	mat[0*2+0] = d = 1.0f / di;
	mat[0*2+1] *= d;
	d = -d;
	mat[1*2+0] *= d;
	d = mat[1*2+0] * di;
	mat[1*2+1] += mat[0*2+1] * d;
	di = mat[1*2+1];
	s *= di;
	mat[1*2+1] = d = 1.0f / di;
	mat[1*2+0] *= d;
	d = -d;
	mat[0*2+1] *= d;
	d = mat[0*2+1] * di;
	mat[0*2+0] += mat[1*2+0] * d;

	return ( s != 0.0f && !FLOAT_IS_NAN( s ) );
#endif
}

/*
=============
CMat2D::toString
=============
*/
const char *CMat2D::toString( int precision ) const {
	return CMyString::floatArraytoString( toFloatPtr(), getDimension(), precision );
}


//===============================================================
//
//	CMat3D
//
//===============================================================

CMat3D mat3_zero( CVec3D( 0, 0, 0 ), CVec3D( 0, 0, 0 ), CVec3D( 0, 0, 0 ) );
CMat3D mat3_identity( CVec3D( 1, 0, 0 ), CVec3D( 0, 1, 0 ), CVec3D( 0, 0, 1 ) );

/*
============
CMat3D::toAngles
============
*/
CEulerAngles CMat3D::toAngles() const {
	CEulerAngles	angles;
	double		theta;
	double		cp;
	float		sp;

	sp = mat[ 0 ][ 2 ];

	// cap off our sin value so that we don't get any NANs
	if ( sp > 1.0f ) {
		sp = 1.0f;
	} else if ( sp < -1.0f ) {
		sp = -1.0f;
	}

	theta = -CMath::asin( sp );
	cp = CMath::cos( theta );

	if ( cp > 8192.0f * CMath::FLT_EPSILON ) {
		angles.pitch	= RAD2DEG( theta );
		angles.yaw		= RAD2DEG( atan2( mat[ 0 ][ 1 ], mat[ 0 ][ 0 ] ) );
		angles.roll		= RAD2DEG( atan2( mat[ 1 ][ 2 ], mat[ 2 ][ 2 ] ) );
	} else {
		angles.pitch	= RAD2DEG( theta );
		angles.yaw		= RAD2DEG( -atan2( mat[ 1 ][ 0 ], mat[ 1 ][ 1 ] ) );
		angles.roll		= 0;
	}
	return angles;
}

/*
============
CMat3D::toQuat
============
*/
CQuaternion CMat3D::toQuat() const {
	CQuaternion		q;
	float		trace;
	float		s;
	float		t;
	int     	i;
	int			j;
	int			k;

	static int 	next[ 3 ] = { 1, 2, 0 };

	trace = mat[ 0 ][ 0 ] + mat[ 1 ][ 1 ] + mat[ 2 ][ 2 ];

	if ( trace > 0.0f ) {

		t = trace + 1.0f;
		s = CMath::invSqrt( t ) * 0.5f;

		q[3] = s * t;
		q[0] = ( mat[ 2 ][ 1 ] - mat[ 1 ][ 2 ] ) * s;
		q[1] = ( mat[ 0 ][ 2 ] - mat[ 2 ][ 0 ] ) * s;
		q[2] = ( mat[ 1 ][ 0 ] - mat[ 0 ][ 1 ] ) * s;

	} else {

		i = 0;
		if ( mat[ 1 ][ 1 ] > mat[ 0 ][ 0 ] ) {
			i = 1;
		}
		if ( mat[ 2 ][ 2 ] > mat[ i ][ i ] ) {
			i = 2;
		}
		j = next[ i ];
		k = next[ j ];

		t = ( mat[ i ][ i ] - ( mat[ j ][ j ] + mat[ k ][ k ] ) ) + 1.0f;
		s = CMath::invSqrt( t ) * 0.5f;

		q[i] = s * t;
		q[3] = ( mat[ k ][ j ] - mat[ j ][ k ] ) * s;
		q[j] = ( mat[ j ][ i ] + mat[ i ][ j ] ) * s;
		q[k] = ( mat[ k ][ i ] + mat[ i ][ k ] ) * s;
	}
	return q;
}

/*
============
CMat3D::toQuat
============
*/
CCompQuaternion CMat3D::toCQuat() const {
	CQuaternion q = toQuat();
	if ( q.w < 0.0f ) {
		return CCompQuaternion( -q.x, -q.y, -q.z );
	}
	return CCompQuaternion( q.x, q.y, q.z );
}

/*
============
CMat3D::toRotation
============
*/
CRotation CMat3D::toRotation() const {
	CRotation	r;
	float		trace;
	float		s;
	float		t;
	int     	i;
	int			j;
	int			k;
	static int 	next[ 3 ] = { 1, 2, 0 };

	trace = mat[ 0 ][ 0 ] + mat[ 1 ][ 1 ] + mat[ 2 ][ 2 ];
	if ( trace > 0.0f ) {

		t = trace + 1.0f;
		s = CMath::invSqrt( t ) * 0.5f;

		r.angle = s * t;
		r.vec[0] = ( mat[ 2 ][ 1 ] - mat[ 1 ][ 2 ] ) * s;
		r.vec[1] = ( mat[ 0 ][ 2 ] - mat[ 2 ][ 0 ] ) * s;
		r.vec[2] = ( mat[ 1 ][ 0 ] - mat[ 0 ][ 1 ] ) * s;

	} else {

		i = 0;
		if ( mat[ 1 ][ 1 ] > mat[ 0 ][ 0 ] ) {
			i = 1;
		}
		if ( mat[ 2 ][ 2 ] > mat[ i ][ i ] ) {
			i = 2;
		}
		j = next[ i ];
		k = next[ j ];

		t = ( mat[ i ][ i ] - ( mat[ j ][ j ] + mat[ k ][ k ] ) ) + 1.0f;
		s = CMath::invSqrt( t ) * 0.5f;

		r.vec[i]	= s * t;
		r.angle		= ( mat[ k ][ j ] - mat[ j ][ k ] ) * s;
		r.vec[j]	= ( mat[ j ][ i ] + mat[ i ][ j ] ) * s;
		r.vec[k]	= ( mat[ k ][ i ] + mat[ i ][ k ] ) * s;
	}
	r.angle = CMath::acos( r.angle );
	if ( CMath::fabs( r.angle ) < 1e-10f ) {
		r.vec.set( 0.0f, 0.0f, 1.0f );
		r.angle = 0.0f;
	} else {
		//vec *= (1.0f / sin( angle ));
		r.vec.toNormal();
		r.vec.fixDegenerateNormal();
		r.angle *= 2.0f * CMath::M_RAD2DEG;
	}

	r.origin.toZero();
	r.axis = *this;
	r.axisValid = true;
	return r;
}

/*
=================
CMat3D::toAngularVelocity
=================
*/
CVec3D CMat3D::toAngularVelocity() const {
	CRotation rotation = toRotation();
	return rotation.getVec() * DEG2RAD( rotation.getAngle() );
}

/*
============
CMat3D::Determinant
============
*/
float CMat3D::determinant() const {

	float det2_12_01 = mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0];
	float det2_12_02 = mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0];
	float det2_12_12 = mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1];

	return mat[0][0] * det2_12_12 - mat[0][1] * det2_12_02 + mat[0][2] * det2_12_01;
}

/*
============
CMat3D::InverseSelf
============
*/
bool CMat3D::inverseSelf() {
	// 18+3+9 = 30 multiplications
	//			 1 division
	CMat3D inverse;
	double det, invDet;

	inverse[0][0] = mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1];
	inverse[1][0] = mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2];
	inverse[2][0] = mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0];

	det = mat[0][0] * inverse[0][0] + mat[0][1] * inverse[1][0] + mat[0][2] * inverse[2][0];

	if ( CMath::fabs( det ) < MATRIX_INVERSE_EPSILON ) {
		return false;
	}

	invDet = 1.0f / det;

	inverse[0][1] = mat[0][2] * mat[2][1] - mat[0][1] * mat[2][2];
	inverse[0][2] = mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1];
	inverse[1][1] = mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0];
	inverse[1][2] = mat[0][2] * mat[1][0] - mat[0][0] * mat[1][2];
	inverse[2][1] = mat[0][1] * mat[2][0] - mat[0][0] * mat[2][1];
	inverse[2][2] = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];

	mat[0][0] = inverse[0][0] * invDet;
	mat[0][1] = inverse[0][1] * invDet;
	mat[0][2] = inverse[0][2] * invDet;

	mat[1][0] = inverse[1][0] * invDet;
	mat[1][1] = inverse[1][1] * invDet;
	mat[1][2] = inverse[1][2] * invDet;

	mat[2][0] = inverse[2][0] * invDet;
	mat[2][1] = inverse[2][1] * invDet;
	mat[2][2] = inverse[2][2] * invDet;

	return true;
}

/*
============
CMat3D::InverseFastSelf
============
*/
bool CMat3D::inverseFastSelf() {
#if 1
	// 18+3+9 = 30 multiplications
	//			 1 division
	CMat3D inverse;
	double det, invDet;

	inverse[0][0] = mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1];
	inverse[1][0] = mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2];
	inverse[2][0] = mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0];

	det = mat[0][0] * inverse[0][0] + mat[0][1] * inverse[1][0] + mat[0][2] * inverse[2][0];

	if ( CMath::fabs( det ) < MATRIX_INVERSE_EPSILON ) {
		return false;
	}

	invDet = 1.0f / det;

	inverse[0][1] = mat[0][2] * mat[2][1] - mat[0][1] * mat[2][2];
	inverse[0][2] = mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1];
	inverse[1][1] = mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0];
	inverse[1][2] = mat[0][2] * mat[1][0] - mat[0][0] * mat[1][2];
	inverse[2][1] = mat[0][1] * mat[2][0] - mat[0][0] * mat[2][1];
	inverse[2][2] = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];

	mat[0][0] = inverse[0][0] * invDet;
	mat[0][1] = inverse[0][1] * invDet;
	mat[0][2] = inverse[0][2] * invDet;

	mat[1][0] = inverse[1][0] * invDet;
	mat[1][1] = inverse[1][1] * invDet;
	mat[1][2] = inverse[1][2] * invDet;

	mat[2][0] = inverse[2][0] * invDet;
	mat[2][1] = inverse[2][1] * invDet;
	mat[2][2] = inverse[2][2] * invDet;

	return true;
#elif 0
	// 3*10 = 30 multiplications
	//		   3 divisions
	float *mat = reinterpret_cast<float *>(this);
	float s;
	double d, di;

	di = mat[0];
	s = di;
	mat[0] = d = 1.0f / di;
	mat[1] *= d;
	mat[2] *= d;
	d = -d;
	mat[3] *= d;
	mat[6] *= d;
	d = mat[3] * di;
	mat[4] += mat[1] * d;
	mat[5] += mat[2] * d;
	d = mat[6] * di;
	mat[7] += mat[1] * d;
	mat[8] += mat[2] * d;
	di = mat[4];
	s *= di;
	mat[4] = d = 1.0f / di;
	mat[3] *= d;
	mat[5] *= d;
	d = -d;
	mat[1] *= d;
	mat[7] *= d;
	d = mat[1] * di;
	mat[0] += mat[3] * d;
	mat[2] += mat[5] * d;
	d = mat[7] * di;
	mat[6] += mat[3] * d;
	mat[8] += mat[5] * d;
	di = mat[8];
	s *= di;
	mat[8] = d = 1.0f / di;
	mat[6] *= d;
	mat[7] *= d;
	d = -d;
	mat[2] *= d;
	mat[5] *= d;
	d = mat[2] * di;
	mat[0] += mat[6] * d;
	mat[1] += mat[7] * d;
	d = mat[5] * di;
	mat[3] += mat[6] * d;
	mat[4] += mat[7] * d;

	return ( s != 0.0f && !FLOAT_IS_NAN( s ) );
#else
	//	4*2+4*4 = 24 multiplications
	//		2*1 =  2 divisions
	CMat2D r0;
	float r1[2], r2[2], r3;
	float det, invDet;
	float *mat = reinterpret_cast<float *>(this);

	// r0 = m0.Inverse();	// 2x2
	det = mat[0*3+0] * mat[1*3+1] - mat[0*3+1] * mat[1*3+0];

	if ( CMath::fabs( det ) < MATRIX_INVERSE_EPSILON ) {
		return false;
	}

	invDet = 1.0f / det;

	r0[0][0] =   mat[1*3+1] * invDet;
	r0[0][1] = - mat[0*3+1] * invDet;
	r0[1][0] = - mat[1*3+0] * invDet;
	r0[1][1] =   mat[0*3+0] * invDet;

	// r1 = r0 * m1;		// 2x1 = 2x2 * 2x1
	r1[0] = r0[0][0] * mat[0*3+2] + r0[0][1] * mat[1*3+2];
	r1[1] = r0[1][0] * mat[0*3+2] + r0[1][1] * mat[1*3+2];

	// r2 = m2 * r1;		// 1x1 = 1x2 * 2x1
	r2[0] = mat[2*3+0] * r1[0] + mat[2*3+1] * r1[1];

	// r3 = r2 - m3;		// 1x1 = 1x1 - 1x1
	r3 = r2[0] - mat[2*3+2];

	// r3.InverseSelf();
	if ( CMath::fabs( r3 ) < MATRIX_INVERSE_EPSILON ) {
		return false;
	}

	r3 = 1.0f / r3;

	// r2 = m2 * r0;		// 1x2 = 1x2 * 2x2
	r2[0] = mat[2*3+0] * r0[0][0] + mat[2*3+1] * r0[1][0];
	r2[1] = mat[2*3+0] * r0[0][1] + mat[2*3+1] * r0[1][1];

	// m2 = r3 * r2;		// 1x2 = 1x1 * 1x2
	mat[2*3+0] = r3 * r2[0];
	mat[2*3+1] = r3 * r2[1];

	// m0 = r0 - r1 * m2;	// 2x2 - 2x1 * 1x2
	mat[0*3+0] = r0[0][0] - r1[0] * mat[2*3+0];
	mat[0*3+1] = r0[0][1] - r1[0] * mat[2*3+1];
	mat[1*3+0] = r0[1][0] - r1[1] * mat[2*3+0];
	mat[1*3+1] = r0[1][1] - r1[1] * mat[2*3+1];

	// m1 = r1 * r3;		// 2x1 = 2x1 * 1x1
	mat[0*3+2] = r1[0] * r3;
	mat[1*3+2] = r1[1] * r3;

	// m3 = -r3;
	mat[2*3+2] = -r3;

	return true;
#endif
}

#if 0
/************************************************************
*
* input:
* mat - pointer to array of 16 floats (source matrix)
* output:
* dst - pointer to array of 16 floats (invert matrix)
*
*************************************************************/
void Invert2(float *mat, float *dst)
{
float tmp[12]; /* temp array for pairs */
float src[16]; /* array of transpose source matrix */
float det; /* determinant */
/* transpose matrix */
for (int i = 0; i < 4; i++) {
src[i] = mat[i*4];
src[i + 4] = mat[i*4 + 1];
src[i + 8] = mat[i*4 + 2];
src[i + 12] = mat[i*4 + 3];
}
/* calculate pairs for first 8 elements (cofactors) */
tmp[0] = src[10] * src[15];
tmp[1] = src[11] * src[14];
tmp[2] = src[9] * src[15];
tmp[3] = src[11] * src[13];
tmp[4] = src[9] * src[14];
tmp[5] = src[10] * src[13];
tmp[6] = src[8] * src[15];
tmp[7] = src[11] * src[12];
tmp[8] = src[8] * src[14];
tmp[9] = src[10] * src[12];
tmp[10] = src[8] * src[13];
tmp[11] = src[9] * src[12];
/* calculate first 8 elements (cofactors) */
dst[0] = tmp[0]*src[5] + tmp[3]*src[6] + tmp[4]*src[7];
dst[0] -= tmp[1]*src[5] + tmp[2]*src[6] + tmp[5]*src[7];
dst[1] = tmp[1]*src[4] + tmp[6]*src[6] + tmp[9]*src[7];
dst[1] -= tmp[0]*src[4] + tmp[7]*src[6] + tmp[8]*src[7];
dst[2] = tmp[2]*src[4] + tmp[7]*src[5] + tmp[10]*src[7];
dst[2] -= tmp[3]*src[4] + tmp[6]*src[5] + tmp[11]*src[7];
dst[3] = tmp[5]*src[4] + tmp[8]*src[5] + tmp[11]*src[6];
dst[3] -= tmp[4]*src[4] + tmp[9]*src[5] + tmp[10]*src[6];
dst[4] = tmp[1]*src[1] + tmp[2]*src[2] + tmp[5]*src[3];
dst[4] -= tmp[0]*src[1] + tmp[3]*src[2] + tmp[4]*src[3];
dst[5] = tmp[0]*src[0] + tmp[7]*src[2] + tmp[8]*src[3];
dst[5] -= tmp[1]*src[0] + tmp[6]*src[2] + tmp[9]*src[3];
dst[6] = tmp[3]*src[0] + tmp[6]*src[1] + tmp[11]*src[3];
dst[6] -= tmp[2]*src[0] + tmp[7]*src[1] + tmp[10]*src[3];
dst[7] = tmp[4]*src[0] + tmp[9]*src[1] + tmp[10]*src[2];
dst[7] -= tmp[5]*src[0] + tmp[8]*src[1] + tmp[11]*src[2];
/* calculate pairs for second 8 elements (cofactors) */
tmp[0] = src[2]*src[7];
tmp[1] = src[3]*src[6];
tmp[2] = src[1]*src[7];
tmp[3] = src[3]*src[5];
tmp[4] = src[1]*src[6];
tmp[5] = src[2]*src[5];
tmp[6] = src[0]*src[7];
tmp[7] = src[3]*src[4];
tmp[8] = src[0]*src[6];
tmp[9] = src[2]*src[4];
tmp[10] = src[0]*src[5];
tmp[11] = src[1]*src[4];
/* calculate second 8 elements (cofactors) */
dst[8] = tmp[0]*src[13] + tmp[3]*src[14] + tmp[4]*src[15];
dst[8] -= tmp[1]*src[13] + tmp[2]*src[14] + tmp[5]*src[15];
dst[9] = tmp[1]*src[12] + tmp[6]*src[14] + tmp[9]*src[15];
dst[9] -= tmp[0]*src[12] + tmp[7]*src[14] + tmp[8]*src[15];
dst[10] = tmp[2]*src[12] + tmp[7]*src[13] + tmp[10]*src[15];
dst[10]-= tmp[3]*src[12] + tmp[6]*src[13] + tmp[11]*src[15];
dst[11] = tmp[5]*src[12] + tmp[8]*src[13] + tmp[11]*src[14];
dst[11]-= tmp[4]*src[12] + tmp[9]*src[13] + tmp[10]*src[14];
dst[12] = tmp[2]*src[10] + tmp[5]*src[11] + tmp[1]*src[9];
dst[12]-= tmp[4]*src[11] + tmp[0]*src[9] + tmp[3]*src[10];
dst[13] = tmp[8]*src[11] + tmp[0]*src[8] + tmp[7]*src[10];
dst[13]-= tmp[6]*src[10] + tmp[9]*src[11] + tmp[1]*src[8];
dst[14] = tmp[6]*src[9] + tmp[11]*src[11] + tmp[3]*src[8];
dst[14]-= tmp[10]*src[11] + tmp[2]*src[8] + tmp[7]*src[9];
dst[15] = tmp[10]*src[10] + tmp[4]*src[8] + tmp[9]*src[9];
dst[15]-= tmp[8]*src[9] + tmp[11]*src[10] + tmp[5]*src[8];
/* calculate determinant */
det=src[0]*dst[0]+src[1]*dst[1]+src[2]*dst[2]+src[3]*dst[3];
/* calculate matrix inverse */
det = 1/det;
for (int j = 0; j < 16; j++)
dst[j] *= det;
}

//==========================================================
/**
Brief description of the program:
1.
Variables (Streaming SIMD Extensions registers) which will contain cofactors and, later, the
lines of the inverted matrix are declared.
2.
Variables which will contain the lines of the reference matrix and, later (after the transposition),
the columns of the original matrix are declared.
3.
Temporary variables and the variable that will contain the matrix determinant are declared.
Streaming SIMD Extensions - Inverse of 4x4 Matrix
9
4 - 11.
Matrix transposition.
12 - 57.
Cofactors calculation. Because in the process of cofactor computation some pairs in threeelement
products are repeated, it is not reasonable to load these pairs anew every time. The
values in the registers with these pairs are formed using shuffle instruction. Cofactors are
calculated row by row (4 elements are placed in 1 SP FP SIMD floating point register).
58 - 63.
Evaluation of determinant and its reciprocal value. 1/det is evaluated using a fast rcpps
command with subsequent approximation using the Newton-Raphson algorithm.
64 - 75.
Multiplication of cofactors by 1/det. Storing the inverse matrix to the address in pointer src.
**/
void PIII_inverse_4x4(float* src)
{
sf_m128 minor0, minor1, minor2, minor3;
sf_m128 row0, row1, row2, row3;
sf_m128 det, tmp1;
tmp1 = _mm_loadh_pi(_mm_loadl_pi(tmp1, (__m64*)(src)), (__m64*)(src+ 4));
row1 = _mm_loadh_pi(_mm_loadl_pi(row1, (__m64*)(src+8)), (__m64*)(src+12));
row0 = _mm_shuffle_ps(tmp1, row1, 0x88);
row1 = _mm_shuffle_ps(row1, tmp1, 0xDD);
tmp1 = _mm_loadh_pi(_mm_loadl_pi(tmp1, (__m64*)(src+ 2)), (__m64*)(src+ 6));
row3 = _mm_loadh_pi(_mm_loadl_pi(row3, (__m64*)(src+10)), (__m64*)(src+14));
row2 = _mm_shuffle_ps(tmp1, row3, 0x88);
row3 = _mm_shuffle_ps(row3, tmp1, 0xDD);
// -----------------------------------------------
tmp1 = _mm_mul_ps(row2, row3);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
minor0 = _mm_mul_ps(row1, tmp1);
minor1 = _mm_mul_ps(row0, tmp1);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor0 = _mm_sub_ps(_mm_mul_ps(row1, tmp1), minor0);
minor1 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), minor1);
minor1 = _mm_shuffle_ps(minor1, minor1, 0x4E);
// -----------------------------------------------
tmp1 = _mm_mul_ps(row1, row2);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
minor0 = _mm_add_ps(_mm_mul_ps(row3, tmp1), minor0);
minor3 = _mm_mul_ps(row0, tmp1);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor0 = _mm_sub_ps(minor0, _mm_mul_ps(row3, tmp1));
minor3 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), minor3);
minor3 = _mm_shuffle_ps(minor3, minor3, 0x4E);
// -----------------------------------------------
tmp1 = _mm_mul_ps(_mm_shuffle_ps(row1, row1, 0x4E), row3);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
row2 = _mm_shuffle_ps(row2, row2, 0x4E);
minor0 = _mm_add_ps(_mm_mul_ps(row2, tmp1), minor0);
minor2 = _mm_mul_ps(row0, tmp1);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor0 = _mm_sub_ps(minor0, _mm_mul_ps(row2, tmp1));
minor2 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), minor2);
minor2 = _mm_shuffle_ps(minor2, minor2, 0x4E);
// -----------------------------------------------
tmp1 = _mm_mul_ps(row0, row1);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
minor2 = _mm_add_ps(_mm_mul_ps(row3, tmp1), minor2);
minor3 = _mm_sub_ps(_mm_mul_ps(row2, tmp1), minor3);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor2 = _mm_sub_ps(_mm_mul_ps(row3, tmp1), minor2);
minor3 = _mm_sub_ps(minor3, _mm_mul_ps(row2, tmp1));
// -----------------------------------------------
tmp1 = _mm_mul_ps(row0, row3);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
minor1 = _mm_sub_ps(minor1, _mm_mul_ps(row2, tmp1));
minor2 = _mm_add_ps(_mm_mul_ps(row1, tmp1), minor2);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor1 = _mm_add_ps(_mm_mul_ps(row2, tmp1), minor1);
minor2 = _mm_sub_ps(minor2, _mm_mul_ps(row1, tmp1));
// -----------------------------------------------
tmp1 = _mm_mul_ps(row0, row2);
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
minor1 = _mm_add_ps(_mm_mul_ps(row3, tmp1), minor1);
minor3 = _mm_sub_ps(minor3, _mm_mul_ps(row1, tmp1));
tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
minor1 = _mm_sub_ps(minor1, _mm_mul_ps(row3, tmp1));
minor3 = _mm_add_ps(_mm_mul_ps(row1, tmp1), minor3);
// -----------------------------------------------
det = _mm_mul_ps(row0, minor0);
det = _mm_add_ps(_mm_shuffle_ps(det, det, 0x4E), det);
det = _mm_add_ss(_mm_shuffle_ps(det, det, 0xB1), det);
tmp1 = _mm_rcp_ss(det);
det = _mm_sub_ss(_mm_add_ss(tmp1, tmp1), _mm_mul_ss(det, _mm_mul_ss(tmp1, tmp1)));
det = _mm_shuffle_ps(det, det, 0x00);
minor0 = _mm_mul_ps(det, minor0);
_mm_storel_pi((__m64*)(src), minor0);
_mm_storeh_pi((__m64*)(src+2), minor0);
minor1 = _mm_mul_ps(det, minor1);
_mm_storel_pi((__m64*)(src+4), minor1);
_mm_storeh_pi((__m64*)(src+6), minor1);
minor2 = _mm_mul_ps(det, minor2);
_mm_storel_pi((__m64*)(src+ 8), minor2);
_mm_storeh_pi((__m64*)(src+10), minor2);
minor3 = _mm_mul_ps(det, minor3);
_mm_storel_pi((__m64*)(src+12), minor3);
_mm_storeh_pi((__m64*)(src+14), minor3);
}
#endif
/*
============
CMat3D::InertiaTranslate
============
*/
CMat3D CMat3D::InertiaTranslate( const float mass, const CVec3D &centerOfMass, const CVec3D &translation ) const {
	CMat3D m;
	CVec3D newCenter;

	newCenter = centerOfMass + translation;

	m[0][0] = mass * ( ( centerOfMass[1] * centerOfMass[1] + centerOfMass[2] * centerOfMass[2] )
				- ( newCenter[1] * newCenter[1] + newCenter[2] * newCenter[2] ) );
	m[1][1] = mass * ( ( centerOfMass[0] * centerOfMass[0] + centerOfMass[2] * centerOfMass[2] )
				- ( newCenter[0] * newCenter[0] + newCenter[2] * newCenter[2] ) );
	m[2][2] = mass * ( ( centerOfMass[0] * centerOfMass[0] + centerOfMass[1] * centerOfMass[1] )
				- ( newCenter[0] * newCenter[0] + newCenter[1] * newCenter[1] ) );

	m[0][1] = m[1][0] = mass * ( newCenter[0] * newCenter[1] - centerOfMass[0] * centerOfMass[1] );
	m[1][2] = m[2][1] = mass * ( newCenter[1] * newCenter[2] - centerOfMass[1] * centerOfMass[2] );
	m[0][2] = m[2][0] = mass * ( newCenter[0] * newCenter[2] - centerOfMass[0] * centerOfMass[2] );

	return (*this) + m;
}

/*
============
CMat3D::InertiaTranslateSelf
============
*/
CMat3D &CMat3D::InertiaTranslateSelf( const float mass, const CVec3D &centerOfMass, const CVec3D &translation ) {
	CMat3D m;
	CVec3D newCenter;

	newCenter = centerOfMass + translation;

	m[0][0] = mass * ( ( centerOfMass[1] * centerOfMass[1] + centerOfMass[2] * centerOfMass[2] )
				- ( newCenter[1] * newCenter[1] + newCenter[2] * newCenter[2] ) );
	m[1][1] = mass * ( ( centerOfMass[0] * centerOfMass[0] + centerOfMass[2] * centerOfMass[2] )
				- ( newCenter[0] * newCenter[0] + newCenter[2] * newCenter[2] ) );
	m[2][2] = mass * ( ( centerOfMass[0] * centerOfMass[0] + centerOfMass[1] * centerOfMass[1] )
				- ( newCenter[0] * newCenter[0] + newCenter[1] * newCenter[1] ) );

	m[0][1] = m[1][0] = mass * ( newCenter[0] * newCenter[1] - centerOfMass[0] * centerOfMass[1] );
	m[1][2] = m[2][1] = mass * ( newCenter[1] * newCenter[2] - centerOfMass[1] * centerOfMass[2] );
	m[0][2] = m[2][0] = mass * ( newCenter[0] * newCenter[2] - centerOfMass[0] * centerOfMass[2] );

	(*this) += m;

	return (*this);
}

/*
============
CMat3D::InertiaRotate
============
*/
CMat3D CMat3D::InertiaRotate( const CMat3D &rotation ) const {
	// NOTE: the rotation matrix is stored column-major
	return rotation.transpose() * (*this) * rotation;
}

/*
============
CMat3D::InertiaRotateSelf
============
*/
CMat3D &CMat3D::InertiaRotateSelf( const CMat3D &rotation ) {
	// NOTE: the rotation matrix is stored column-major
	*this = rotation.transpose() * (*this) * rotation;
	return *this;
}



void CMat3D::batchTransform(CVec3D *pointArray, int numPoints) const
{
	SMF_ASSERT(pointArray || numPoints == 0);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (!pointArray)
		return;
#endif
	for(int i = 0; i < numPoints; ++i)
		pointArray[i] = *this * pointArray[i];
}

void CMat3D::batchTransform(CVec3D *pointArray, int numPoints, int stride) const
{
	SMF_ASSERT(pointArray || numPoints == 0);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (!pointArray)
		return;
#endif
	SMF_ASSERT(stride >= (int)sizeof(CVec3D));
	sf_u8 *data = reinterpret_cast<sf_u8*>(pointArray);
	for(int i = 0; i < numPoints; ++i)
	{
		CVec3D *v = reinterpret_cast<CVec3D*>(data + stride*i);
		*v = *this * *v;
	}
}

void CMat3D::batchTransform(CVec4D *vectorArray, int numVectors) const
{
	SMF_ASSERT(vectorArray || numVectors == 0);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (!vectorArray)
		return;
#endif
	for(int i = 0; i < numVectors; ++i)
		vectorArray[i] = *this * vectorArray[i];
}

void CMat3D::batchTransform(CVec4D *vectorArray, int numVectors, int stride) const
{
	SMF_ASSERT(vectorArray || numVectors == 0);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (!vectorArray)
		return;
#endif
	SMF_ASSERT(stride >= (int)sizeof(CVec4D));
	sf_u8 *data = reinterpret_cast<sf_u8*>(vectorArray);
	for(int i = 0; i < numVectors; ++i)
	{
		CVec4D *v = reinterpret_cast<CVec4D*>(data + stride*i);
		*v = (*this) * (*v);
	}
}



/*
=============
CMat3D::toString
=============
*/
const char *CMat3D::toString( int precision ) const {
	return CMyString::floatArraytoString( toFloatPtr(), getDimension(), precision );
}

CMat3D CMat3D::operator *(const CQuaternion &quat) const
{
	return (*this) * quat.toMat3();
}

CMat4D CMat4D::operator *(const CQuaternion &rhs) const
{

	CMat3D rot=rhs.toMat3();
	return *this * rot;

}

CVec3D CMat3D::transform(const CVec3D &vector) const
{
	return transform(vector.x, vector.y, vector.z);
}

CVec3D CMat3D::transformLeft(const CVec3D &vector) const
{
	return CVec3D((vector*Col(0)),
	              (vector*Col(1)),
	              (vector*Col(2)));
}

CVec3D CMat3D::transform(float x, float y, float z) const
{
	return CVec3D(Row(0)*CVec3D( x,y,z),
	             Row(1)*CVec3D( x,y,z),
	             Row(2)*CVec3D( x,y,z));
}

CVec4D CMat3D::transform(const CVec4D &vector) const
{
	return CVec4D((Row(0)* vector.toVec3()),
	              (Row(1)* vector.toVec3()),
	              (Row(2)* vector.toVec3()),
	              vector.w);
}

bool CMat3D::solveAxb(CVec3D b, CVec3D &x) const
{
	// solve by pivotization.
	float v00 = mat[0][0];
	float v10 = mat[0][1];
	float v20 = mat[0][2];

	float v01 = mat[1][0];
	float v11 = mat[1][1];
	float v21 = mat[1][2];

	float v02 = mat[2][0];
	float v12 = mat[2][1];
	float v22 = mat[2][2];

	float av00 = CMath::fabs(v00);
	float av10 = CMath::fabs(v10);
	float av20 = CMath::fabs(v20);

	// find which item in first column has largest absolute value.
	if (av10 >= av00 && av10 >= av20)
	{
		Swap(v00, v10);
		Swap(v01, v11);
		Swap(v02, v12);
		Swap(b[0], b[1]);
	}
	else if (v20 >= v00)
	{
		Swap(v00, v20);
		Swap(v01, v21);
		Swap(v02, v22);
		Swap(b[0], b[2]);
	}

	/* a b c | x
	   d e f | y
	   g h i | z , where |a| >= |d| && |a| >= |g| */

	if (CMath::equalsAbs(v00, 0.f))
		return false;

	// scale row so that leading element is one.
	float denom = 1.f / v00;
//	v00 = 1.f;
	v01 *= denom;
	v02 *= denom;
	b[0] *= denom;

	/* 1 b c | x
	   d e f | y
	   g h i | z */

	// zero first column of second and third rows.
	v11 -= v10 * v01;
	v12 -= v10 * v02;
	b[1] -= v10 * b[0];

	v21 -= v20 * v01;
	v22 -= v20 * v02;
	b[2] -= v20 * b[0];

	/* 1 b c | x
	   0 e f | y
	   0 h i | z */

	// Pivotize again.
	if (CMath::fabs(v21) >CMath:: abs(v11))
	{
		Swap(v11, v21);
		Swap(v12, v22);
		Swap(b[1], b[2]);
	}

	if (CMath::equalsAbs(v11, 0.f))
		return false;

	/* 1 b c | x
	   0 e f | y
	   0 h i | z, where |e| >= |h| */

	denom = 1.f / v11;
//	v11 = 1.f;
	v12 *= denom;
	b[1] *= denom;

	/* 1 b c | x
	   0 1 f | y
	   0 h i | z */

	v22 -= v21 * v12;
	b[2] -= v21 * b[1];

	/* 1 b c | x
	   0 1 f | y
	   0 0 i | z */

	if (CMath::equalsAbs(v22, 0.f))
		return false;

	x[2] = b[2] / v22;
	x[1] = b[1] - x[2] * v12;
	x[0] = b[0] - x[2] * v02 - x[1] * v01;

	return true;
}

void CMat3D::scaleRow(int row, float scalar)
{
	SMF_ASSERT(row >= 0);
	SMF_ASSERT(row < 3);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (row < 0 || row >= 3)
		return;
#endif
	SMF_ASSERT(MATH::isFinite(scalar));
	Row(row) *= scalar;
}

void CMat3D::scaleCol(int col, float scalar)
{
	SMF_ASSERT(col >= 0);
	SMF_ASSERT(col < 3);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (col < 0 || col >= 3)
		return;
#endif
	mat[col][0] *= scalar;
	mat[col][1] *= scalar;
	mat[col][2] *= scalar;
}


//===============================================================
//
//	CMat4D
//
//===============================================================



CMat4D mat4_zero( CVec4D( 0, 0, 0, 0 ), CVec4D( 0, 0, 0, 0 ), CVec4D( 0, 0, 0, 0 ), CVec4D( 0, 0, 0, 0 ) );
CMat4D mat4_identity( CVec4D( 1, 0, 0, 0 ), CVec4D( 0, 1, 0, 0 ), CVec4D( 0, 0, 1, 0 ), CVec4D( 0, 0, 0, 1 ) );

CMat4D::CMat4D(const CQuaternion  &orientation)
{
if (CMath::MATH_AUTOMATIC_SSE) {
	CSIMD::getProcessor()->quat_to_mat4x4(orientation.q, sseSignMask4, Row(0).toM128Ptr());
}else{
	CSIMD::getGenProcessor()->quat_to_mat4x4(orientation.q, sseSignMask4,Row(0).toM128Ptr());

//	setRotatePart(orientation);
//	setRow(3, 0, 0, 0, 1);
//	setCol3(3, 0, 0, 0);

}
}

CMat4D::CMat4D(const CQuaternion  &orientation, const CVec3D &translation)
{
if (CMath::MATH_AUTOMATIC_SSE) {
	CSIMD::getProcessor()->quat_to_mat4x4(orientation.q, CVec4D(translation, 1.f).v, Row(0).toM128Ptr());
}else{
	CSIMD::getGenProcessor()->quat_to_mat4x4(orientation.q, CVec4D(translation, 1.f).v, Row(0).toM128Ptr());

}
}

CMat4D CMat4D::operator *(const CMatJoint3x4 &rhs) const
{
	///\todo SSE.
	CMat4D r;
	CVec3D col0=rhs.Col(0);
	CVec3D col1=rhs.Col(1);
	CVec3D col2=rhs.Col(2);
	CVec3D col3=rhs.Col(3);

	r.mat[0][0] = mat[0].toVec3() * col0; //DOT3STRIDED(v[0], c0, 4);
	r.mat[0][1] = mat[0].toVec3() * col1; //DOT3STRIDED(v[0], c1, 4);
	r.mat[0][2] = mat[0].toVec3() * col2; //DOT3STRIDED(v[0], c2, 4);
	r.mat[0][3] = (mat[0].toVec3() * col3)+mat[0][3]; //DOT3STRIDED(v[0], c3, 4) + v[0][3];

	r.mat[1][0] = mat[1].toVec3() * col0; // DOT3STRIDED(v[1], c0, 4);
	r.mat[1][1] = mat[1].toVec3() * col1; // DOT3STRIDED(v[1], c1, 4);
	r.mat[1][2] = mat[1].toVec3() * col2; // DOT3STRIDED(v[1], c2, 4);
	r.mat[1][3] = (mat[1].toVec3() * col3)+mat[1][3]; // DOT3STRIDED(v[1], c3, 4) + v[1][3];

	r.mat[2][0] = mat[2].toVec3() * col0; // DOT3STRIDED(v[2], c0, 4);
	r.mat[2][1] = mat[2].toVec3() * col1; // DOT3STRIDED(v[2], c1, 4);
	r.mat[2][2] = mat[2].toVec3() * col2; // DOT3STRIDED(v[2], c2, 4);
	r.mat[2][3] = (mat[2].toVec3() * col3)+mat[2][3]; // DOT3STRIDED(v[2], c3, 4) + v[2][3];

	r.mat[3][0] = mat[3].toVec3() * col0; // DOT3STRIDED(v[3], c0, 4);
	r.mat[3][1] = mat[3].toVec3() * col1; // DOT3STRIDED(v[3], c1, 4);
	r.mat[3][2] = mat[3].toVec3() * col2; // DOT3STRIDED(v[3], c2, 4);
	r.mat[3][3] = (mat[3].toVec3() * col3)+mat[3][3]; // DOT3STRIDED(v[3], c3, 4) + v[3][3];

	return r;
}


/*
============
CMat4D::transpose
============
*/
CMat4D CMat4D::transpose() const {
	CMat4D	transpose;
	int		i, j;

	for( i = 0; i < 4; i++ ) {
		for( j = 0; j < 4; j++ ) {
			transpose[ i ][ j ] = mat[ j ][ i ];
        }
	}
	return transpose;
}

/*
============
CMat4D::transposeSelf
============
*/
CMat4D &CMat4D::transposeSelf() {
	float	temp;
	int		i, j;

	for( i = 0; i < 4; i++ ) {
		for( j = i + 1; j < 4; j++ ) {
			temp = mat[ i ][ j ];
			mat[ i ][ j ] = mat[ j ][ i ];
			mat[ j ][ i ] = temp;
        }
	}
	return *this;
}

/*
============
CMat4D::Determinant
============
*/
float CMat4D::determinant() const {
	if(CMath::MATH_AUTOMATIC_SSE) {
		return CSIMD::getProcessor()->mat4x4_determinant(this);
	}else{
		return CSIMD::getGenProcessor()->mat4x4_determinant(this);

#if 0
	// 2x2 sub-determinants
	float det2_01_01 = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
	float det2_01_02 = mat[0][0] * mat[1][2] - mat[0][2] * mat[1][0];
	float det2_01_03 = mat[0][0] * mat[1][3] - mat[0][3] * mat[1][0];
	float det2_01_12 = mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1];
	float det2_01_13 = mat[0][1] * mat[1][3] - mat[0][3] * mat[1][1];
	float det2_01_23 = mat[0][2] * mat[1][3] - mat[0][3] * mat[1][2];

	// 3x3 sub-determinants
	float det3_201_012 = mat[2][0] * det2_01_12 - mat[2][1] * det2_01_02 + mat[2][2] * det2_01_01;
	float det3_201_013 = mat[2][0] * det2_01_13 - mat[2][1] * det2_01_03 + mat[2][3] * det2_01_01;
	float det3_201_023 = mat[2][0] * det2_01_23 - mat[2][2] * det2_01_03 + mat[2][3] * det2_01_02;
	float det3_201_123 = mat[2][1] * det2_01_23 - mat[2][2] * det2_01_13 + mat[2][3] * det2_01_12;

	return ( - det3_201_123 * mat[3][0] + det3_201_023 * mat[3][1] - det3_201_013 * mat[3][2] + det3_201_012 * mat[3][3] );
#endif
	}
}

CMat4D CMat4D::inverse() const throw(Exception::CMathException){


	CMat4D invMat;

	invMat = *this;
	int r = invMat.inverseSelf();
	SMF_ASSERT(r !=0)


	return invMat;
}

/*
============
CMat4D::InverseSelf
============
*/
bool CMat4D::inverseSelf() {

	if(CMath::MATH_AUTOMATIC_SSE){
		CSIMD::getProcessor()->mat4x4_inverse(this, this);
		return true;
	}else{
		CSIMD::getGenProcessor()->mat4x4_inverse(this, this);
		return true;
	}
}

/*
============
CMat4D::InverseFastSelf
============
*/
bool CMat4D::inverseFastSelf() {
#if 0 //def MATX_SIMD
	SIMDProcessor->matX_inverse_4x4(this->toFloatPtr());

#elif 0
	// 84+4+16 = 104 multiplications
	//			   1 division
	double det, invDet;

	// 2x2 sub-determinants required to calculate 4x4 determinant
	float det2_01_01 = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
	float det2_01_02 = mat[0][0] * mat[1][2] - mat[0][2] * mat[1][0];
	float det2_01_03 = mat[0][0] * mat[1][3] - mat[0][3] * mat[1][0];
	float det2_01_12 = mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1];
	float det2_01_13 = mat[0][1] * mat[1][3] - mat[0][3] * mat[1][1];
	float det2_01_23 = mat[0][2] * mat[1][3] - mat[0][3] * mat[1][2];

	// 3x3 sub-determinants required to calculate 4x4 determinant
	float det3_201_012 = mat[2][0] * det2_01_12 - mat[2][1] * det2_01_02 + mat[2][2] * det2_01_01;
	float det3_201_013 = mat[2][0] * det2_01_13 - mat[2][1] * det2_01_03 + mat[2][3] * det2_01_01;
	float det3_201_023 = mat[2][0] * det2_01_23 - mat[2][2] * det2_01_03 + mat[2][3] * det2_01_02;
	float det3_201_123 = mat[2][1] * det2_01_23 - mat[2][2] * det2_01_13 + mat[2][3] * det2_01_12;

	det = ( - det3_201_123 * mat[3][0] + det3_201_023 * mat[3][1] - det3_201_013 * mat[3][2] + det3_201_012 * mat[3][3] );

	if ( CMath::fabs( det ) < MATRIX_INVERSE_EPSILON ) {
		return false;
	}

	invDet = 1.0f / det;

	// remaining 2x2 sub-determinants
	float det2_03_01 = mat[0][0] * mat[3][1] - mat[0][1] * mat[3][0];
	float det2_03_02 = mat[0][0] * mat[3][2] - mat[0][2] * mat[3][0];
	float det2_03_03 = mat[0][0] * mat[3][3] - mat[0][3] * mat[3][0];
	float det2_03_12 = mat[0][1] * mat[3][2] - mat[0][2] * mat[3][1];
	float det2_03_13 = mat[0][1] * mat[3][3] - mat[0][3] * mat[3][1];
	float det2_03_23 = mat[0][2] * mat[3][3] - mat[0][3] * mat[3][2];

	float det2_13_01 = mat[1][0] * mat[3][1] - mat[1][1] * mat[3][0];
	float det2_13_02 = mat[1][0] * mat[3][2] - mat[1][2] * mat[3][0];
	float det2_13_03 = mat[1][0] * mat[3][3] - mat[1][3] * mat[3][0];
	float det2_13_12 = mat[1][1] * mat[3][2] - mat[1][2] * mat[3][1];
	float det2_13_13 = mat[1][1] * mat[3][3] - mat[1][3] * mat[3][1];
	float det2_13_23 = mat[1][2] * mat[3][3] - mat[1][3] * mat[3][2];

	// remaining 3x3 sub-determinants
	float det3_203_012 = mat[2][0] * det2_03_12 - mat[2][1] * det2_03_02 + mat[2][2] * det2_03_01;
	float det3_203_013 = mat[2][0] * det2_03_13 - mat[2][1] * det2_03_03 + mat[2][3] * det2_03_01;
	float det3_203_023 = mat[2][0] * det2_03_23 - mat[2][2] * det2_03_03 + mat[2][3] * det2_03_02;
	float det3_203_123 = mat[2][1] * det2_03_23 - mat[2][2] * det2_03_13 + mat[2][3] * det2_03_12;

	float det3_213_012 = mat[2][0] * det2_13_12 - mat[2][1] * det2_13_02 + mat[2][2] * det2_13_01;
	float det3_213_013 = mat[2][0] * det2_13_13 - mat[2][1] * det2_13_03 + mat[2][3] * det2_13_01;
	float det3_213_023 = mat[2][0] * det2_13_23 - mat[2][2] * det2_13_03 + mat[2][3] * det2_13_02;
	float det3_213_123 = mat[2][1] * det2_13_23 - mat[2][2] * det2_13_13 + mat[2][3] * det2_13_12;

	float det3_301_012 = mat[3][0] * det2_01_12 - mat[3][1] * det2_01_02 + mat[3][2] * det2_01_01;
	float det3_301_013 = mat[3][0] * det2_01_13 - mat[3][1] * det2_01_03 + mat[3][3] * det2_01_01;
	float det3_301_023 = mat[3][0] * det2_01_23 - mat[3][2] * det2_01_03 + mat[3][3] * det2_01_02;
	float det3_301_123 = mat[3][1] * det2_01_23 - mat[3][2] * det2_01_13 + mat[3][3] * det2_01_12;

	mat[0][0] =	- det3_213_123 * invDet;
	mat[1][0] = + det3_213_023 * invDet;
	mat[2][0] = - det3_213_013 * invDet;
	mat[3][0] = + det3_213_012 * invDet;

	mat[0][1] = + det3_203_123 * invDet;
	mat[1][1] = - det3_203_023 * invDet;
	mat[2][1] = + det3_203_013 * invDet;
	mat[3][1] = - det3_203_012 * invDet;

	mat[0][2] = + det3_301_123 * invDet;
	mat[1][2] = - det3_301_023 * invDet;
	mat[2][2] = + det3_301_013 * invDet;
	mat[3][2] = - det3_301_012 * invDet;

	mat[0][3] = - det3_201_123 * invDet;
	mat[1][3] = + det3_201_023 * invDet;
	mat[2][3] = - det3_201_013 * invDet;
	mat[3][3] = + det3_201_012 * invDet;

	return true;
#elif 0
	// 4*18 = 72 multiplications
	//		   4 divisions
	float *mat = reinterpret_cast<float *>(this);
	float s;
	double d, di;

	di = mat[0];
	s = di;
	mat[0] = d = 1.0f / di;
	mat[1] *= d;
	mat[2] *= d;
	mat[3] *= d;
	d = -d;
	mat[4] *= d;
	mat[8] *= d;
	mat[12] *= d;
	d = mat[4] * di;
	mat[5] += mat[1] * d;
	mat[6] += mat[2] * d;
	mat[7] += mat[3] * d;
	d = mat[8] * di;
	mat[9] += mat[1] * d;
	mat[10] += mat[2] * d;
	mat[11] += mat[3] * d;
	d = mat[12] * di;
	mat[13] += mat[1] * d;
	mat[14] += mat[2] * d;
	mat[15] += mat[3] * d;
	di = mat[5];
	s *= di;
	mat[5] = d = 1.0f / di;
	mat[4] *= d;
	mat[6] *= d;
	mat[7] *= d;
	d = -d;
	mat[1] *= d;
	mat[9] *= d;
	mat[13] *= d;
	d = mat[1] * di;
	mat[0] += mat[4] * d;
	mat[2] += mat[6] * d;
	mat[3] += mat[7] * d;
	d = mat[9] * di;
	mat[8] += mat[4] * d;
	mat[10] += mat[6] * d;
	mat[11] += mat[7] * d;
	d = mat[13] * di;
	mat[12] += mat[4] * d;
	mat[14] += mat[6] * d;
	mat[15] += mat[7] * d;
	di = mat[10];
	s *= di;
	mat[10] = d = 1.0f / di;
	mat[8] *= d;
	mat[9] *= d;
	mat[11] *= d;
	d = -d;
	mat[2] *= d;
	mat[6] *= d;
	mat[14] *= d;
	d = mat[2] * di;
	mat[0] += mat[8] * d;
	mat[1] += mat[9] * d;
	mat[3] += mat[11] * d;
	d = mat[6] * di;
	mat[4] += mat[8] * d;
	mat[5] += mat[9] * d;
	mat[7] += mat[11] * d;
	d = mat[14] * di;
	mat[12] += mat[8] * d;
	mat[13] += mat[9] * d;
	mat[15] += mat[11] * d;
	di = mat[15];
	s *= di;
	mat[15] = d = 1.0f / di;
	mat[12] *= d;
	mat[13] *= d;
	mat[14] *= d;
	d = -d;
	mat[3] *= d;
	mat[7] *= d;
	mat[11] *= d;
	d = mat[3] * di;
	mat[0] += mat[12] * d;
	mat[1] += mat[13] * d;
	mat[2] += mat[14] * d;
	d = mat[7] * di;
	mat[4] += mat[12] * d;
	mat[5] += mat[13] * d;
	mat[6] += mat[14] * d;
	d = mat[11] * di;
	mat[8] += mat[12] * d;
	mat[9] += mat[13] * d;
	mat[10] += mat[14] * d;

	return ( s != 0.0f && !FLOAT_IS_NAN( s ) );
#else
	//	6*8+2*6 = 60 multiplications
	//		2*1 =  2 divisions
	CMat2D r0, r1, r2, r3;
	float a, det, invDet;
	float *mat = reinterpret_cast<float *>(this);

	// r0 = m0.Inverse();
	det = mat[0*4+0] * mat[1*4+1] - mat[0*4+1] * mat[1*4+0];

	if ( CMath::fabs( det ) < MATRIX_INVERSE_EPSILON ) {
		return false;
	}

	invDet = 1.0f / det;

	r0[0][0] =   mat[1*4+1] * invDet;
	r0[0][1] = - mat[0*4+1] * invDet;
	r0[1][0] = - mat[1*4+0] * invDet;
	r0[1][1] =   mat[0*4+0] * invDet;

	// r1 = r0 * m1;
	r1[0][0] = r0[0][0] * mat[0*4+2] + r0[0][1] * mat[1*4+2];
	r1[0][1] = r0[0][0] * mat[0*4+3] + r0[0][1] * mat[1*4+3];
	r1[1][0] = r0[1][0] * mat[0*4+2] + r0[1][1] * mat[1*4+2];
	r1[1][1] = r0[1][0] * mat[0*4+3] + r0[1][1] * mat[1*4+3];

	// r2 = m2 * r1;
	r2[0][0] = mat[2*4+0] * r1[0][0] + mat[2*4+1] * r1[1][0];
	r2[0][1] = mat[2*4+0] * r1[0][1] + mat[2*4+1] * r1[1][1];
	r2[1][0] = mat[3*4+0] * r1[0][0] + mat[3*4+1] * r1[1][0];
	r2[1][1] = mat[3*4+0] * r1[0][1] + mat[3*4+1] * r1[1][1];

	// r3 = r2 - m3;
	r3[0][0] = r2[0][0] - mat[2*4+2];
	r3[0][1] = r2[0][1] - mat[2*4+3];
	r3[1][0] = r2[1][0] - mat[3*4+2];
	r3[1][1] = r2[1][1] - mat[3*4+3];

	// r3.InverseSelf();
	det = r3[0][0] * r3[1][1] - r3[0][1] * r3[1][0];

	if ( CMath::fabs( det ) < MATRIX_INVERSE_EPSILON ) {
		return false;
	}

	invDet = 1.0f / det;

	a = r3[0][0];
	r3[0][0] =   r3[1][1] * invDet;
	r3[0][1] = - r3[0][1] * invDet;
	r3[1][0] = - r3[1][0] * invDet;
	r3[1][1] =   a * invDet;

	// r2 = m2 * r0;
	r2[0][0] = mat[2*4+0] * r0[0][0] + mat[2*4+1] * r0[1][0];
	r2[0][1] = mat[2*4+0] * r0[0][1] + mat[2*4+1] * r0[1][1];
	r2[1][0] = mat[3*4+0] * r0[0][0] + mat[3*4+1] * r0[1][0];
	r2[1][1] = mat[3*4+0] * r0[0][1] + mat[3*4+1] * r0[1][1];

	// m2 = r3 * r2;
	mat[2*4+0] = r3[0][0] * r2[0][0] + r3[0][1] * r2[1][0];
	mat[2*4+1] = r3[0][0] * r2[0][1] + r3[0][1] * r2[1][1];
	mat[3*4+0] = r3[1][0] * r2[0][0] + r3[1][1] * r2[1][0];
	mat[3*4+1] = r3[1][0] * r2[0][1] + r3[1][1] * r2[1][1];

	// m0 = r0 - r1 * m2;
	mat[0*4+0] = r0[0][0] - r1[0][0] * mat[2*4+0] - r1[0][1] * mat[3*4+0];
	mat[0*4+1] = r0[0][1] - r1[0][0] * mat[2*4+1] - r1[0][1] * mat[3*4+1];
	mat[1*4+0] = r0[1][0] - r1[1][0] * mat[2*4+0] - r1[1][1] * mat[3*4+0];
	mat[1*4+1] = r0[1][1] - r1[1][0] * mat[2*4+1] - r1[1][1] * mat[3*4+1];

	// m1 = r1 * r3;
	mat[0*4+2] = r1[0][0] * r3[0][0] + r1[0][1] * r3[1][0];
	mat[0*4+3] = r1[0][0] * r3[0][1] + r1[0][1] * r3[1][1];
	mat[1*4+2] = r1[1][0] * r3[0][0] + r1[1][1] * r3[1][0];
	mat[1*4+3] = r1[1][0] * r3[0][1] + r1[1][1] * r3[1][1];

	// m3 = -r3;
	mat[2*4+2] = -r3[0][0];
	mat[2*4+3] = -r3[0][1];
	mat[3*4+2] = -r3[1][0];
	mat[3*4+3] = -r3[1][1];

	return true;
#endif
}

float CMat4D::Determinant3() const
{
	SMF_ASSERT(CMat3DPart().isFinite());
if (CMath::MATH_AUTOMATIC_SSE) {
	return CSIMD::getProcessor()->mat3x4_determinant(Row(0).toM128Ptr());
}else{
	return CSIMD::getGenProcessor()->mat3x4_determinant(Row(0).toM128Ptr());
	}
}
CVec3D CMat4D::extractScale() const
{
	return CVec3D(Col3(0).getLenght(), Col3(1).getLenght(), Col3(2).getLenght());
}




const char *CMat4D::toString( int precision ) const {
	return CMyString::floatArraytoString( toFloatPtr(), getDimension(), precision );
}


void CMat4D::set3x3Part(const CMat3D &r)
{
	SMF_ASSERT(r.isFinite());
	mat[0][0] = r[0][0]; mat[0][1] = r[0][1]; mat[0][2] = r[0][2];
	mat[1][0] = r[1][0]; mat[1][1] = r[1][1]; mat[1][2] = r[1][2];
	mat[2][0] = r[2][0]; mat[2][1] = r[2][1]; mat[2][2] = r[2][2];
}

void CMat4D::set3x4Part(const CMatJoint3x4 &r)
{
	SMF_ASSERT(r.isFinite());

		mat[0][0] = r[0*4+0]; mat[0][1] = r[0*4+1]; mat[0][2] = r[0*4+2]; mat[0][3] = r[0*4+3];
		mat[1][0] = r[1*4+0]; mat[1][1] = r[1*4+1]; mat[1][2] = r[1*4+2]; mat[1][3] = r[1*4+3];
		mat[2][0] = r[2*4+0]; mat[2][1] = r[2*4+1]; mat[2][2] = r[2*4+2]; mat[2][3] = r[2*4+3];

}
void CMat4D::setCol3(int column, const CVec3D &columnVector)
{
	setCol3(column, columnVector.x, columnVector.y, columnVector.z);
}

void CMat4D::setCol3(int column, const float *data)
{
	SMF_ASSERT(data);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (!data)
		return;
#endif
	setCol3(column, data[0], data[1], data[2]);
}

void CMat4D::setCol3(int column, float m_0c, float m_1c, float m_2c)
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
	mat[0][column] = m_0c;
	mat[1][column] = m_1c;
	mat[2][column] = m_2c;
}

void CMat4D::set(float _00, float _01, float _02, float _03,
				   float _10, float _11, float _12, float _13,
				   float _20, float _21, float _22, float _23,
				   float _30, float _31, float _32, float _33)
{
if(CMath::MATH_AUTOMATIC_SSE){
	mat[0] = MATH::SET_PS(_03, _02, _01, _00);
	mat[1] = MATH::SET_PS(_13, _12, _11, _10);
	mat[2] = MATH::SET_PS(_23, _22, _21, _20);
	mat[3] = MATH::SET_PS(_33, _32, _31, _30);

}else{
	mat[0][0] = _00; mat[0][1] = _01; mat[0][2] = _02; mat[0][3] = _03;
	mat[1][0] = _10; mat[1][1] = _11; mat[1][2] = _12; mat[1][3] = _13;
	mat[2][0] = _20; mat[2][1] = _21; mat[2][2] = _22; mat[2][3] = _23;
	mat[3][0] = _30; mat[3][1] = _31; mat[3][2] = _32; mat[3][3] = _33;
}
}


void CMat4D::setTranslatePart(float tx, float ty, float tz)
{
	setCol3(3, tx, ty, tz);
}

void CMat4D::setTranslatePart(const CVec3D &offset)
{
	setCol3(3, offset);
}

void CMat4D::setRotatePartX(float angle)
{
	Set3x3PartRotateX(*this, angle);
}

void CMat4D::setRotatePartY(float angle)
{
	Set3x3PartRotateY(*this, angle);
}

void CMat4D::setRotatePartZ(float angle)
{
	Set3x3PartRotateZ(*this, angle);
}

void CMat4D::setRotatePart(const CVec3D &a, float angle)
{
	SMF_ASSERT(a.isNormalized());
	SMF_ASSERT(MATH::isFinite(angle));

	const float c = CMath::cos(angle);
	const float c1 = (1.f-c);
	const float s = CMath::sin(angle);

	mat[0][0] = c+c1*a.x*a.x;
	mat[1][0] = c1*a.x*a.y+s*a.z;
	mat[2][0] = c1*a.x*a.z-s*a.y;

	mat[0][1] = c1*a.x*a.y-s*a.z;
	mat[1][1] = c+c1*a.y*a.y;
	mat[2][1] = c1*a.y*a.z+s*a.x;

	mat[0][2] = c1*a.x*a.z+s*a.y;
	mat[1][2] = c1*a.y*a.z-s*a.x;
	mat[2][2] = c+c1*a.z*a.z;
}

void CMat4D::setRotatePart(const CQuaternion &q)
{
	SetMatrixRotatePart(*this, q);
}

void CMat4D::setRow(int row, const CVec3D &rowVector, float m_r3)
{
	setRow(row, rowVector.x, rowVector.y, rowVector.z, m_r3);
}

void CMat4D::setRow(int row, const CVec4D &rowVector)
{
	setRow(row, rowVector.x, rowVector.y, rowVector.z, rowVector.w);
}

void CMat4D::setRow(int row, const float *data)
{
	SMF_ASSERT(data);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (!data)
		return;
#endif
	setRow(row, data[0], data[1], data[2], data[3]);
}

void CMat4D::setRow(int row, float m_r0, float m_r1, float m_r2, float m_r3)
{
	SMF_ASSERT(row >= 0);
	SMF_ASSERT(row < 4);
	SMF_ASSERT(MATH::isFinite(m_r0));
	SMF_ASSERT(MATH::isFinite(m_r1));
	SMF_ASSERT(MATH::isFinite(m_r2));
	SMF_ASSERT(MATH::isFinite(m_r3));
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (row < 0 || row >= 4)
		return; // Benign failure
#endif

if (CMath::MATH_AUTOMATIC_SSE) {
	mat[row].v = MATH::SET_PS(m_r3, m_r2, m_r1, m_r0);
}else{
	mat[row][0] = m_r0;
	mat[row][1] = m_r1;
	mat[row][2] = m_r2;
	mat[row][3] = m_r3;
}
}

bool CMat4D::hasUnitaryScale(float epsilon) const
{
	CVec3D scale = extractScale();
	return scale.compare(CVec3D(1.f, 1.f, 1.f), epsilon);
}

bool CMat4D::hasNegativeScale() const
{
	return Determinant3() < 0.f;
}

bool CMat4D::hasUniformScale(float epsilon) const
{
	CVec3D scale = extractScale();
	return CMath::equalsAbs(scale.x, scale.y, epsilon) && CMath::equalsAbs(scale.x, scale.z, epsilon);
}

bool CMat4D::isRowOrthogonal3(float epsilon) const
{
	return Row3(0).isPerpendicular(Row3(1), epsilon)
		&& Row3(0).isPerpendicular(Row3(2), epsilon)
		&& Row3(1).isPerpendicular(Row3(2), epsilon);
}

bool CMat4D::isColOrthogonal3(float epsilon) const
{
	return Col3(0).isPerpendicular(Col3(1), epsilon)
		&& Col3(0).isPerpendicular(Col3(2), epsilon)
		&& Col3(1).isPerpendicular(Col3(2), epsilon);
}

bool CMat4D::IsOrthonormal3(float epsilon) const
{
	///\todo Epsilon magnitudes don't match.
	return isColOrthogonal3(epsilon) && Row3(0).isNormalized(epsilon) && Row3(1).isNormalized(epsilon) && Row3(2).isNormalized(epsilon);
}


//===============================================================
//
//	CMat5D
//
//===============================================================

CMat5D mat5_zero( CVec5D( 0, 0, 0, 0, 0 ), CVec5D( 0, 0, 0, 0, 0 ), CVec5D( 0, 0, 0, 0, 0 ), CVec5D( 0, 0, 0, 0, 0 ), CVec5D( 0, 0, 0, 0, 0 ) );
CMat5D mat5_identity( CVec5D( 1, 0, 0, 0, 0 ), CVec5D( 0, 1, 0, 0, 0 ), CVec5D( 0, 0, 1, 0, 0 ), CVec5D( 0, 0, 0, 1, 0 ), CVec5D( 0, 0, 0, 0, 1 ) );

/*
============
CMat5D::transpose
============
*/
CMat5D CMat5D::transpose() const {
	CMat5D	transpose;
	int		i, j;

	for( i = 0; i < 5; i++ ) {
		for( j = 0; j < 5; j++ ) {
			transpose[ i ][ j ] = mat[ j ][ i ];
        }
	}
	return transpose;
}

/*
============
CMat5D::transposeSelf
============
*/
CMat5D &CMat5D::transposeSelf() {
	float	temp;
	int		i, j;

	for( i = 0; i < 5; i++ ) {
		for( j = i + 1; j < 5; j++ ) {
			temp = mat[ i ][ j ];
			mat[ i ][ j ] = mat[ j ][ i ];
			mat[ j ][ i ] = temp;
        }
	}
	return *this;
}

/*
============
CMat5D::Determinant
============
*/
float CMat5D::determinant() const {

	// 2x2 sub-determinants required to calculate 5x5 determinant
	float det2_34_01 = mat[3][0] * mat[4][1] - mat[3][1] * mat[4][0];
	float det2_34_02 = mat[3][0] * mat[4][2] - mat[3][2] * mat[4][0];
	float det2_34_03 = mat[3][0] * mat[4][3] - mat[3][3] * mat[4][0];
	float det2_34_04 = mat[3][0] * mat[4][4] - mat[3][4] * mat[4][0];
	float det2_34_12 = mat[3][1] * mat[4][2] - mat[3][2] * mat[4][1];
	float det2_34_13 = mat[3][1] * mat[4][3] - mat[3][3] * mat[4][1];
	float det2_34_14 = mat[3][1] * mat[4][4] - mat[3][4] * mat[4][1];
	float det2_34_23 = mat[3][2] * mat[4][3] - mat[3][3] * mat[4][2];
	float det2_34_24 = mat[3][2] * mat[4][4] - mat[3][4] * mat[4][2];
	float det2_34_34 = mat[3][3] * mat[4][4] - mat[3][4] * mat[4][3];

	// 3x3 sub-determinants required to calculate 5x5 determinant
	float det3_234_012 = mat[2][0] * det2_34_12 - mat[2][1] * det2_34_02 + mat[2][2] * det2_34_01;
	float det3_234_013 = mat[2][0] * det2_34_13 - mat[2][1] * det2_34_03 + mat[2][3] * det2_34_01;
	float det3_234_014 = mat[2][0] * det2_34_14 - mat[2][1] * det2_34_04 + mat[2][4] * det2_34_01;
	float det3_234_023 = mat[2][0] * det2_34_23 - mat[2][2] * det2_34_03 + mat[2][3] * det2_34_02;
	float det3_234_024 = mat[2][0] * det2_34_24 - mat[2][2] * det2_34_04 + mat[2][4] * det2_34_02;
	float det3_234_034 = mat[2][0] * det2_34_34 - mat[2][3] * det2_34_04 + mat[2][4] * det2_34_03;
	float det3_234_123 = mat[2][1] * det2_34_23 - mat[2][2] * det2_34_13 + mat[2][3] * det2_34_12;
	float det3_234_124 = mat[2][1] * det2_34_24 - mat[2][2] * det2_34_14 + mat[2][4] * det2_34_12;
	float det3_234_134 = mat[2][1] * det2_34_34 - mat[2][3] * det2_34_14 + mat[2][4] * det2_34_13;
	float det3_234_234 = mat[2][2] * det2_34_34 - mat[2][3] * det2_34_24 + mat[2][4] * det2_34_23;

	// 4x4 sub-determinants required to calculate 5x5 determinant
	float det4_1234_0123 = mat[1][0] * det3_234_123 - mat[1][1] * det3_234_023 + mat[1][2] * det3_234_013 - mat[1][3] * det3_234_012;
	float det4_1234_0124 = mat[1][0] * det3_234_124 - mat[1][1] * det3_234_024 + mat[1][2] * det3_234_014 - mat[1][4] * det3_234_012;
	float det4_1234_0134 = mat[1][0] * det3_234_134 - mat[1][1] * det3_234_034 + mat[1][3] * det3_234_014 - mat[1][4] * det3_234_013;
	float det4_1234_0234 = mat[1][0] * det3_234_234 - mat[1][2] * det3_234_034 + mat[1][3] * det3_234_024 - mat[1][4] * det3_234_023;
	float det4_1234_1234 = mat[1][1] * det3_234_234 - mat[1][2] * det3_234_134 + mat[1][3] * det3_234_124 - mat[1][4] * det3_234_123;

	// determinant of 5x5 matrix
	return mat[0][0] * det4_1234_1234 - mat[0][1] * det4_1234_0234 + mat[0][2] * det4_1234_0134 - mat[0][3] * det4_1234_0124 + mat[0][4] * det4_1234_0123;
}

/*
============
CMat5D::InverseSelf
============
*/
bool CMat5D::inverseSelf() {
	// 280+5+25 = 310 multiplications
	//				1 division
	double det, invDet;

	// 2x2 sub-determinants required to calculate 5x5 determinant
	float det2_34_01 = mat[3][0] * mat[4][1] - mat[3][1] * mat[4][0];
	float det2_34_02 = mat[3][0] * mat[4][2] - mat[3][2] * mat[4][0];
	float det2_34_03 = mat[3][0] * mat[4][3] - mat[3][3] * mat[4][0];
	float det2_34_04 = mat[3][0] * mat[4][4] - mat[3][4] * mat[4][0];
	float det2_34_12 = mat[3][1] * mat[4][2] - mat[3][2] * mat[4][1];
	float det2_34_13 = mat[3][1] * mat[4][3] - mat[3][3] * mat[4][1];
	float det2_34_14 = mat[3][1] * mat[4][4] - mat[3][4] * mat[4][1];
	float det2_34_23 = mat[3][2] * mat[4][3] - mat[3][3] * mat[4][2];
	float det2_34_24 = mat[3][2] * mat[4][4] - mat[3][4] * mat[4][2];
	float det2_34_34 = mat[3][3] * mat[4][4] - mat[3][4] * mat[4][3];

	// 3x3 sub-determinants required to calculate 5x5 determinant
	float det3_234_012 = mat[2][0] * det2_34_12 - mat[2][1] * det2_34_02 + mat[2][2] * det2_34_01;
	float det3_234_013 = mat[2][0] * det2_34_13 - mat[2][1] * det2_34_03 + mat[2][3] * det2_34_01;
	float det3_234_014 = mat[2][0] * det2_34_14 - mat[2][1] * det2_34_04 + mat[2][4] * det2_34_01;
	float det3_234_023 = mat[2][0] * det2_34_23 - mat[2][2] * det2_34_03 + mat[2][3] * det2_34_02;
	float det3_234_024 = mat[2][0] * det2_34_24 - mat[2][2] * det2_34_04 + mat[2][4] * det2_34_02;
	float det3_234_034 = mat[2][0] * det2_34_34 - mat[2][3] * det2_34_04 + mat[2][4] * det2_34_03;
	float det3_234_123 = mat[2][1] * det2_34_23 - mat[2][2] * det2_34_13 + mat[2][3] * det2_34_12;
	float det3_234_124 = mat[2][1] * det2_34_24 - mat[2][2] * det2_34_14 + mat[2][4] * det2_34_12;
	float det3_234_134 = mat[2][1] * det2_34_34 - mat[2][3] * det2_34_14 + mat[2][4] * det2_34_13;
	float det3_234_234 = mat[2][2] * det2_34_34 - mat[2][3] * det2_34_24 + mat[2][4] * det2_34_23;

	// 4x4 sub-determinants required to calculate 5x5 determinant
	float det4_1234_0123 = mat[1][0] * det3_234_123 - mat[1][1] * det3_234_023 + mat[1][2] * det3_234_013 - mat[1][3] * det3_234_012;
	float det4_1234_0124 = mat[1][0] * det3_234_124 - mat[1][1] * det3_234_024 + mat[1][2] * det3_234_014 - mat[1][4] * det3_234_012;
	float det4_1234_0134 = mat[1][0] * det3_234_134 - mat[1][1] * det3_234_034 + mat[1][3] * det3_234_014 - mat[1][4] * det3_234_013;
	float det4_1234_0234 = mat[1][0] * det3_234_234 - mat[1][2] * det3_234_034 + mat[1][3] * det3_234_024 - mat[1][4] * det3_234_023;
	float det4_1234_1234 = mat[1][1] * det3_234_234 - mat[1][2] * det3_234_134 + mat[1][3] * det3_234_124 - mat[1][4] * det3_234_123;

	// determinant of 5x5 matrix
	det = mat[0][0] * det4_1234_1234 - mat[0][1] * det4_1234_0234 + mat[0][2] * det4_1234_0134 - mat[0][3] * det4_1234_0124 + mat[0][4] * det4_1234_0123;

	if( CMath::fabs( det ) < MATRIX_INVERSE_EPSILON ) {
		return false;
	}

	invDet = 1.0f / det;

	// remaining 2x2 sub-determinants
	float det2_23_01 = mat[2][0] * mat[3][1] - mat[2][1] * mat[3][0];
	float det2_23_02 = mat[2][0] * mat[3][2] - mat[2][2] * mat[3][0];
	float det2_23_03 = mat[2][0] * mat[3][3] - mat[2][3] * mat[3][0];
	float det2_23_04 = mat[2][0] * mat[3][4] - mat[2][4] * mat[3][0];
	float det2_23_12 = mat[2][1] * mat[3][2] - mat[2][2] * mat[3][1];
	float det2_23_13 = mat[2][1] * mat[3][3] - mat[2][3] * mat[3][1];
	float det2_23_14 = mat[2][1] * mat[3][4] - mat[2][4] * mat[3][1];
	float det2_23_23 = mat[2][2] * mat[3][3] - mat[2][3] * mat[3][2];
	float det2_23_24 = mat[2][2] * mat[3][4] - mat[2][4] * mat[3][2];
	float det2_23_34 = mat[2][3] * mat[3][4] - mat[2][4] * mat[3][3];
	float det2_24_01 = mat[2][0] * mat[4][1] - mat[2][1] * mat[4][0];
	float det2_24_02 = mat[2][0] * mat[4][2] - mat[2][2] * mat[4][0];
	float det2_24_03 = mat[2][0] * mat[4][3] - mat[2][3] * mat[4][0];
	float det2_24_04 = mat[2][0] * mat[4][4] - mat[2][4] * mat[4][0];
	float det2_24_12 = mat[2][1] * mat[4][2] - mat[2][2] * mat[4][1];
	float det2_24_13 = mat[2][1] * mat[4][3] - mat[2][3] * mat[4][1];
	float det2_24_14 = mat[2][1] * mat[4][4] - mat[2][4] * mat[4][1];
	float det2_24_23 = mat[2][2] * mat[4][3] - mat[2][3] * mat[4][2];
	float det2_24_24 = mat[2][2] * mat[4][4] - mat[2][4] * mat[4][2];
	float det2_24_34 = mat[2][3] * mat[4][4] - mat[2][4] * mat[4][3];

	// remaining 3x3 sub-determinants
	float det3_123_012 = mat[1][0] * det2_23_12 - mat[1][1] * det2_23_02 + mat[1][2] * det2_23_01;
	float det3_123_013 = mat[1][0] * det2_23_13 - mat[1][1] * det2_23_03 + mat[1][3] * det2_23_01;
	float det3_123_014 = mat[1][0] * det2_23_14 - mat[1][1] * det2_23_04 + mat[1][4] * det2_23_01;
	float det3_123_023 = mat[1][0] * det2_23_23 - mat[1][2] * det2_23_03 + mat[1][3] * det2_23_02;
	float det3_123_024 = mat[1][0] * det2_23_24 - mat[1][2] * det2_23_04 + mat[1][4] * det2_23_02;
	float det3_123_034 = mat[1][0] * det2_23_34 - mat[1][3] * det2_23_04 + mat[1][4] * det2_23_03;
	float det3_123_123 = mat[1][1] * det2_23_23 - mat[1][2] * det2_23_13 + mat[1][3] * det2_23_12;
	float det3_123_124 = mat[1][1] * det2_23_24 - mat[1][2] * det2_23_14 + mat[1][4] * det2_23_12;
	float det3_123_134 = mat[1][1] * det2_23_34 - mat[1][3] * det2_23_14 + mat[1][4] * det2_23_13;
	float det3_123_234 = mat[1][2] * det2_23_34 - mat[1][3] * det2_23_24 + mat[1][4] * det2_23_23;
	float det3_124_012 = mat[1][0] * det2_24_12 - mat[1][1] * det2_24_02 + mat[1][2] * det2_24_01;
	float det3_124_013 = mat[1][0] * det2_24_13 - mat[1][1] * det2_24_03 + mat[1][3] * det2_24_01;
	float det3_124_014 = mat[1][0] * det2_24_14 - mat[1][1] * det2_24_04 + mat[1][4] * det2_24_01;
	float det3_124_023 = mat[1][0] * det2_24_23 - mat[1][2] * det2_24_03 + mat[1][3] * det2_24_02;
	float det3_124_024 = mat[1][0] * det2_24_24 - mat[1][2] * det2_24_04 + mat[1][4] * det2_24_02;
	float det3_124_034 = mat[1][0] * det2_24_34 - mat[1][3] * det2_24_04 + mat[1][4] * det2_24_03;
	float det3_124_123 = mat[1][1] * det2_24_23 - mat[1][2] * det2_24_13 + mat[1][3] * det2_24_12;
	float det3_124_124 = mat[1][1] * det2_24_24 - mat[1][2] * det2_24_14 + mat[1][4] * det2_24_12;
	float det3_124_134 = mat[1][1] * det2_24_34 - mat[1][3] * det2_24_14 + mat[1][4] * det2_24_13;
	float det3_124_234 = mat[1][2] * det2_24_34 - mat[1][3] * det2_24_24 + mat[1][4] * det2_24_23;
	float det3_134_012 = mat[1][0] * det2_34_12 - mat[1][1] * det2_34_02 + mat[1][2] * det2_34_01;
	float det3_134_013 = mat[1][0] * det2_34_13 - mat[1][1] * det2_34_03 + mat[1][3] * det2_34_01;
	float det3_134_014 = mat[1][0] * det2_34_14 - mat[1][1] * det2_34_04 + mat[1][4] * det2_34_01;
	float det3_134_023 = mat[1][0] * det2_34_23 - mat[1][2] * det2_34_03 + mat[1][3] * det2_34_02;
	float det3_134_024 = mat[1][0] * det2_34_24 - mat[1][2] * det2_34_04 + mat[1][4] * det2_34_02;
	float det3_134_034 = mat[1][0] * det2_34_34 - mat[1][3] * det2_34_04 + mat[1][4] * det2_34_03;
	float det3_134_123 = mat[1][1] * det2_34_23 - mat[1][2] * det2_34_13 + mat[1][3] * det2_34_12;
	float det3_134_124 = mat[1][1] * det2_34_24 - mat[1][2] * det2_34_14 + mat[1][4] * det2_34_12;
	float det3_134_134 = mat[1][1] * det2_34_34 - mat[1][3] * det2_34_14 + mat[1][4] * det2_34_13;
	float det3_134_234 = mat[1][2] * det2_34_34 - mat[1][3] * det2_34_24 + mat[1][4] * det2_34_23;

	// remaining 4x4 sub-determinants
	float det4_0123_0123 = mat[0][0] * det3_123_123 - mat[0][1] * det3_123_023 + mat[0][2] * det3_123_013 - mat[0][3] * det3_123_012;
	float det4_0123_0124 = mat[0][0] * det3_123_124 - mat[0][1] * det3_123_024 + mat[0][2] * det3_123_014 - mat[0][4] * det3_123_012;
	float det4_0123_0134 = mat[0][0] * det3_123_134 - mat[0][1] * det3_123_034 + mat[0][3] * det3_123_014 - mat[0][4] * det3_123_013;
	float det4_0123_0234 = mat[0][0] * det3_123_234 - mat[0][2] * det3_123_034 + mat[0][3] * det3_123_024 - mat[0][4] * det3_123_023;
	float det4_0123_1234 = mat[0][1] * det3_123_234 - mat[0][2] * det3_123_134 + mat[0][3] * det3_123_124 - mat[0][4] * det3_123_123;
	float det4_0124_0123 = mat[0][0] * det3_124_123 - mat[0][1] * det3_124_023 + mat[0][2] * det3_124_013 - mat[0][3] * det3_124_012;
	float det4_0124_0124 = mat[0][0] * det3_124_124 - mat[0][1] * det3_124_024 + mat[0][2] * det3_124_014 - mat[0][4] * det3_124_012;
	float det4_0124_0134 = mat[0][0] * det3_124_134 - mat[0][1] * det3_124_034 + mat[0][3] * det3_124_014 - mat[0][4] * det3_124_013;
	float det4_0124_0234 = mat[0][0] * det3_124_234 - mat[0][2] * det3_124_034 + mat[0][3] * det3_124_024 - mat[0][4] * det3_124_023;
	float det4_0124_1234 = mat[0][1] * det3_124_234 - mat[0][2] * det3_124_134 + mat[0][3] * det3_124_124 - mat[0][4] * det3_124_123;
	float det4_0134_0123 = mat[0][0] * det3_134_123 - mat[0][1] * det3_134_023 + mat[0][2] * det3_134_013 - mat[0][3] * det3_134_012;
	float det4_0134_0124 = mat[0][0] * det3_134_124 - mat[0][1] * det3_134_024 + mat[0][2] * det3_134_014 - mat[0][4] * det3_134_012;
	float det4_0134_0134 = mat[0][0] * det3_134_134 - mat[0][1] * det3_134_034 + mat[0][3] * det3_134_014 - mat[0][4] * det3_134_013;
	float det4_0134_0234 = mat[0][0] * det3_134_234 - mat[0][2] * det3_134_034 + mat[0][3] * det3_134_024 - mat[0][4] * det3_134_023;
	float det4_0134_1234 = mat[0][1] * det3_134_234 - mat[0][2] * det3_134_134 + mat[0][3] * det3_134_124 - mat[0][4] * det3_134_123;
	float det4_0234_0123 = mat[0][0] * det3_234_123 - mat[0][1] * det3_234_023 + mat[0][2] * det3_234_013 - mat[0][3] * det3_234_012;
	float det4_0234_0124 = mat[0][0] * det3_234_124 - mat[0][1] * det3_234_024 + mat[0][2] * det3_234_014 - mat[0][4] * det3_234_012;
	float det4_0234_0134 = mat[0][0] * det3_234_134 - mat[0][1] * det3_234_034 + mat[0][3] * det3_234_014 - mat[0][4] * det3_234_013;
	float det4_0234_0234 = mat[0][0] * det3_234_234 - mat[0][2] * det3_234_034 + mat[0][3] * det3_234_024 - mat[0][4] * det3_234_023;
	float det4_0234_1234 = mat[0][1] * det3_234_234 - mat[0][2] * det3_234_134 + mat[0][3] * det3_234_124 - mat[0][4] * det3_234_123;

	mat[0][0] =  det4_1234_1234 * invDet;
	mat[0][1] = -det4_0234_1234 * invDet;
	mat[0][2] =  det4_0134_1234 * invDet;
	mat[0][3] = -det4_0124_1234 * invDet;
	mat[0][4] =  det4_0123_1234 * invDet;

	mat[1][0] = -det4_1234_0234 * invDet;
	mat[1][1] =  det4_0234_0234 * invDet;
	mat[1][2] = -det4_0134_0234 * invDet;
	mat[1][3] =  det4_0124_0234 * invDet;
	mat[1][4] = -det4_0123_0234 * invDet;

	mat[2][0] =  det4_1234_0134 * invDet;
	mat[2][1] = -det4_0234_0134 * invDet;
	mat[2][2] =  det4_0134_0134 * invDet;
	mat[2][3] = -det4_0124_0134 * invDet;
	mat[2][4] =  det4_0123_0134 * invDet;

	mat[3][0] = -det4_1234_0124 * invDet;
	mat[3][1] =  det4_0234_0124 * invDet;
	mat[3][2] = -det4_0134_0124 * invDet;
	mat[3][3] =  det4_0124_0124 * invDet;
	mat[3][4] = -det4_0123_0124 * invDet;

	mat[4][0] =  det4_1234_0123 * invDet;
	mat[4][1] = -det4_0234_0123 * invDet;
	mat[4][2] =  det4_0134_0123 * invDet;
	mat[4][3] = -det4_0124_0123 * invDet;
	mat[4][4] =  det4_0123_0123 * invDet;

	return true;
}

/*
============
CMat5D::InverseFastSelf
============
*/
bool CMat5D::inverseFastSelf() {
#if 0
	// 280+5+25 = 310 multiplications
	//				1 division
	double det, invDet;

	// 2x2 sub-determinants required to calculate 5x5 determinant
	float det2_34_01 = mat[3][0] * mat[4][1] - mat[3][1] * mat[4][0];
	float det2_34_02 = mat[3][0] * mat[4][2] - mat[3][2] * mat[4][0];
	float det2_34_03 = mat[3][0] * mat[4][3] - mat[3][3] * mat[4][0];
	float det2_34_04 = mat[3][0] * mat[4][4] - mat[3][4] * mat[4][0];
	float det2_34_12 = mat[3][1] * mat[4][2] - mat[3][2] * mat[4][1];
	float det2_34_13 = mat[3][1] * mat[4][3] - mat[3][3] * mat[4][1];
	float det2_34_14 = mat[3][1] * mat[4][4] - mat[3][4] * mat[4][1];
	float det2_34_23 = mat[3][2] * mat[4][3] - mat[3][3] * mat[4][2];
	float det2_34_24 = mat[3][2] * mat[4][4] - mat[3][4] * mat[4][2];
	float det2_34_34 = mat[3][3] * mat[4][4] - mat[3][4] * mat[4][3];

	// 3x3 sub-determinants required to calculate 5x5 determinant
	float det3_234_012 = mat[2][0] * det2_34_12 - mat[2][1] * det2_34_02 + mat[2][2] * det2_34_01;
	float det3_234_013 = mat[2][0] * det2_34_13 - mat[2][1] * det2_34_03 + mat[2][3] * det2_34_01;
	float det3_234_014 = mat[2][0] * det2_34_14 - mat[2][1] * det2_34_04 + mat[2][4] * det2_34_01;
	float det3_234_023 = mat[2][0] * det2_34_23 - mat[2][2] * det2_34_03 + mat[2][3] * det2_34_02;
	float det3_234_024 = mat[2][0] * det2_34_24 - mat[2][2] * det2_34_04 + mat[2][4] * det2_34_02;
	float det3_234_034 = mat[2][0] * det2_34_34 - mat[2][3] * det2_34_04 + mat[2][4] * det2_34_03;
	float det3_234_123 = mat[2][1] * det2_34_23 - mat[2][2] * det2_34_13 + mat[2][3] * det2_34_12;
	float det3_234_124 = mat[2][1] * det2_34_24 - mat[2][2] * det2_34_14 + mat[2][4] * det2_34_12;
	float det3_234_134 = mat[2][1] * det2_34_34 - mat[2][3] * det2_34_14 + mat[2][4] * det2_34_13;
	float det3_234_234 = mat[2][2] * det2_34_34 - mat[2][3] * det2_34_24 + mat[2][4] * det2_34_23;

	// 4x4 sub-determinants required to calculate 5x5 determinant
	float det4_1234_0123 = mat[1][0] * det3_234_123 - mat[1][1] * det3_234_023 + mat[1][2] * det3_234_013 - mat[1][3] * det3_234_012;
	float det4_1234_0124 = mat[1][0] * det3_234_124 - mat[1][1] * det3_234_024 + mat[1][2] * det3_234_014 - mat[1][4] * det3_234_012;
	float det4_1234_0134 = mat[1][0] * det3_234_134 - mat[1][1] * det3_234_034 + mat[1][3] * det3_234_014 - mat[1][4] * det3_234_013;
	float det4_1234_0234 = mat[1][0] * det3_234_234 - mat[1][2] * det3_234_034 + mat[1][3] * det3_234_024 - mat[1][4] * det3_234_023;
	float det4_1234_1234 = mat[1][1] * det3_234_234 - mat[1][2] * det3_234_134 + mat[1][3] * det3_234_124 - mat[1][4] * det3_234_123;

	// determinant of 5x5 matrix
	det = mat[0][0] * det4_1234_1234 - mat[0][1] * det4_1234_0234 + mat[0][2] * det4_1234_0134 - mat[0][3] * det4_1234_0124 + mat[0][4] * det4_1234_0123;

	if( CMath::fabs( det ) < MATRIX_INVERSE_EPSILON ) {
		return false;
	}

	invDet = 1.0f / det;

	// remaining 2x2 sub-determinants
	float det2_23_01 = mat[2][0] * mat[3][1] - mat[2][1] * mat[3][0];
	float det2_23_02 = mat[2][0] * mat[3][2] - mat[2][2] * mat[3][0];
	float det2_23_03 = mat[2][0] * mat[3][3] - mat[2][3] * mat[3][0];
	float det2_23_04 = mat[2][0] * mat[3][4] - mat[2][4] * mat[3][0];
	float det2_23_12 = mat[2][1] * mat[3][2] - mat[2][2] * mat[3][1];
	float det2_23_13 = mat[2][1] * mat[3][3] - mat[2][3] * mat[3][1];
	float det2_23_14 = mat[2][1] * mat[3][4] - mat[2][4] * mat[3][1];
	float det2_23_23 = mat[2][2] * mat[3][3] - mat[2][3] * mat[3][2];
	float det2_23_24 = mat[2][2] * mat[3][4] - mat[2][4] * mat[3][2];
	float det2_23_34 = mat[2][3] * mat[3][4] - mat[2][4] * mat[3][3];
	float det2_24_01 = mat[2][0] * mat[4][1] - mat[2][1] * mat[4][0];
	float det2_24_02 = mat[2][0] * mat[4][2] - mat[2][2] * mat[4][0];
	float det2_24_03 = mat[2][0] * mat[4][3] - mat[2][3] * mat[4][0];
	float det2_24_04 = mat[2][0] * mat[4][4] - mat[2][4] * mat[4][0];
	float det2_24_12 = mat[2][1] * mat[4][2] - mat[2][2] * mat[4][1];
	float det2_24_13 = mat[2][1] * mat[4][3] - mat[2][3] * mat[4][1];
	float det2_24_14 = mat[2][1] * mat[4][4] - mat[2][4] * mat[4][1];
	float det2_24_23 = mat[2][2] * mat[4][3] - mat[2][3] * mat[4][2];
	float det2_24_24 = mat[2][2] * mat[4][4] - mat[2][4] * mat[4][2];
	float det2_24_34 = mat[2][3] * mat[4][4] - mat[2][4] * mat[4][3];

	// remaining 3x3 sub-determinants
	float det3_123_012 = mat[1][0] * det2_23_12 - mat[1][1] * det2_23_02 + mat[1][2] * det2_23_01;
	float det3_123_013 = mat[1][0] * det2_23_13 - mat[1][1] * det2_23_03 + mat[1][3] * det2_23_01;
	float det3_123_014 = mat[1][0] * det2_23_14 - mat[1][1] * det2_23_04 + mat[1][4] * det2_23_01;
	float det3_123_023 = mat[1][0] * det2_23_23 - mat[1][2] * det2_23_03 + mat[1][3] * det2_23_02;
	float det3_123_024 = mat[1][0] * det2_23_24 - mat[1][2] * det2_23_04 + mat[1][4] * det2_23_02;
	float det3_123_034 = mat[1][0] * det2_23_34 - mat[1][3] * det2_23_04 + mat[1][4] * det2_23_03;
	float det3_123_123 = mat[1][1] * det2_23_23 - mat[1][2] * det2_23_13 + mat[1][3] * det2_23_12;
	float det3_123_124 = mat[1][1] * det2_23_24 - mat[1][2] * det2_23_14 + mat[1][4] * det2_23_12;
	float det3_123_134 = mat[1][1] * det2_23_34 - mat[1][3] * det2_23_14 + mat[1][4] * det2_23_13;
	float det3_123_234 = mat[1][2] * det2_23_34 - mat[1][3] * det2_23_24 + mat[1][4] * det2_23_23;
	float det3_124_012 = mat[1][0] * det2_24_12 - mat[1][1] * det2_24_02 + mat[1][2] * det2_24_01;
	float det3_124_013 = mat[1][0] * det2_24_13 - mat[1][1] * det2_24_03 + mat[1][3] * det2_24_01;
	float det3_124_014 = mat[1][0] * det2_24_14 - mat[1][1] * det2_24_04 + mat[1][4] * det2_24_01;
	float det3_124_023 = mat[1][0] * det2_24_23 - mat[1][2] * det2_24_03 + mat[1][3] * det2_24_02;
	float det3_124_024 = mat[1][0] * det2_24_24 - mat[1][2] * det2_24_04 + mat[1][4] * det2_24_02;
	float det3_124_034 = mat[1][0] * det2_24_34 - mat[1][3] * det2_24_04 + mat[1][4] * det2_24_03;
	float det3_124_123 = mat[1][1] * det2_24_23 - mat[1][2] * det2_24_13 + mat[1][3] * det2_24_12;
	float det3_124_124 = mat[1][1] * det2_24_24 - mat[1][2] * det2_24_14 + mat[1][4] * det2_24_12;
	float det3_124_134 = mat[1][1] * det2_24_34 - mat[1][3] * det2_24_14 + mat[1][4] * det2_24_13;
	float det3_124_234 = mat[1][2] * det2_24_34 - mat[1][3] * det2_24_24 + mat[1][4] * det2_24_23;
	float det3_134_012 = mat[1][0] * det2_34_12 - mat[1][1] * det2_34_02 + mat[1][2] * det2_34_01;
	float det3_134_013 = mat[1][0] * det2_34_13 - mat[1][1] * det2_34_03 + mat[1][3] * det2_34_01;
	float det3_134_014 = mat[1][0] * det2_34_14 - mat[1][1] * det2_34_04 + mat[1][4] * det2_34_01;
	float det3_134_023 = mat[1][0] * det2_34_23 - mat[1][2] * det2_34_03 + mat[1][3] * det2_34_02;
	float det3_134_024 = mat[1][0] * det2_34_24 - mat[1][2] * det2_34_04 + mat[1][4] * det2_34_02;
	float det3_134_034 = mat[1][0] * det2_34_34 - mat[1][3] * det2_34_04 + mat[1][4] * det2_34_03;
	float det3_134_123 = mat[1][1] * det2_34_23 - mat[1][2] * det2_34_13 + mat[1][3] * det2_34_12;
	float det3_134_124 = mat[1][1] * det2_34_24 - mat[1][2] * det2_34_14 + mat[1][4] * det2_34_12;
	float det3_134_134 = mat[1][1] * det2_34_34 - mat[1][3] * det2_34_14 + mat[1][4] * det2_34_13;
	float det3_134_234 = mat[1][2] * det2_34_34 - mat[1][3] * det2_34_24 + mat[1][4] * det2_34_23;

	// remaining 4x4 sub-determinants
	float det4_0123_0123 = mat[0][0] * det3_123_123 - mat[0][1] * det3_123_023 + mat[0][2] * det3_123_013 - mat[0][3] * det3_123_012;
	float det4_0123_0124 = mat[0][0] * det3_123_124 - mat[0][1] * det3_123_024 + mat[0][2] * det3_123_014 - mat[0][4] * det3_123_012;
	float det4_0123_0134 = mat[0][0] * det3_123_134 - mat[0][1] * det3_123_034 + mat[0][3] * det3_123_014 - mat[0][4] * det3_123_013;
	float det4_0123_0234 = mat[0][0] * det3_123_234 - mat[0][2] * det3_123_034 + mat[0][3] * det3_123_024 - mat[0][4] * det3_123_023;
	float det4_0123_1234 = mat[0][1] * det3_123_234 - mat[0][2] * det3_123_134 + mat[0][3] * det3_123_124 - mat[0][4] * det3_123_123;
	float det4_0124_0123 = mat[0][0] * det3_124_123 - mat[0][1] * det3_124_023 + mat[0][2] * det3_124_013 - mat[0][3] * det3_124_012;
	float det4_0124_0124 = mat[0][0] * det3_124_124 - mat[0][1] * det3_124_024 + mat[0][2] * det3_124_014 - mat[0][4] * det3_124_012;
	float det4_0124_0134 = mat[0][0] * det3_124_134 - mat[0][1] * det3_124_034 + mat[0][3] * det3_124_014 - mat[0][4] * det3_124_013;
	float det4_0124_0234 = mat[0][0] * det3_124_234 - mat[0][2] * det3_124_034 + mat[0][3] * det3_124_024 - mat[0][4] * det3_124_023;
	float det4_0124_1234 = mat[0][1] * det3_124_234 - mat[0][2] * det3_124_134 + mat[0][3] * det3_124_124 - mat[0][4] * det3_124_123;
	float det4_0134_0123 = mat[0][0] * det3_134_123 - mat[0][1] * det3_134_023 + mat[0][2] * det3_134_013 - mat[0][3] * det3_134_012;
	float det4_0134_0124 = mat[0][0] * det3_134_124 - mat[0][1] * det3_134_024 + mat[0][2] * det3_134_014 - mat[0][4] * det3_134_012;
	float det4_0134_0134 = mat[0][0] * det3_134_134 - mat[0][1] * det3_134_034 + mat[0][3] * det3_134_014 - mat[0][4] * det3_134_013;
	float det4_0134_0234 = mat[0][0] * det3_134_234 - mat[0][2] * det3_134_034 + mat[0][3] * det3_134_024 - mat[0][4] * det3_134_023;
	float det4_0134_1234 = mat[0][1] * det3_134_234 - mat[0][2] * det3_134_134 + mat[0][3] * det3_134_124 - mat[0][4] * det3_134_123;
	float det4_0234_0123 = mat[0][0] * det3_234_123 - mat[0][1] * det3_234_023 + mat[0][2] * det3_234_013 - mat[0][3] * det3_234_012;
	float det4_0234_0124 = mat[0][0] * det3_234_124 - mat[0][1] * det3_234_024 + mat[0][2] * det3_234_014 - mat[0][4] * det3_234_012;
	float det4_0234_0134 = mat[0][0] * det3_234_134 - mat[0][1] * det3_234_034 + mat[0][3] * det3_234_014 - mat[0][4] * det3_234_013;
	float det4_0234_0234 = mat[0][0] * det3_234_234 - mat[0][2] * det3_234_034 + mat[0][3] * det3_234_024 - mat[0][4] * det3_234_023;
	float det4_0234_1234 = mat[0][1] * det3_234_234 - mat[0][2] * det3_234_134 + mat[0][3] * det3_234_124 - mat[0][4] * det3_234_123;

	mat[0][0] =  det4_1234_1234 * invDet;
	mat[0][1] = -det4_0234_1234 * invDet;
	mat[0][2] =  det4_0134_1234 * invDet;
	mat[0][3] = -det4_0124_1234 * invDet;
	mat[0][4] =  det4_0123_1234 * invDet;

	mat[1][0] = -det4_1234_0234 * invDet;
	mat[1][1] =  det4_0234_0234 * invDet;
	mat[1][2] = -det4_0134_0234 * invDet;
	mat[1][3] =  det4_0124_0234 * invDet;
	mat[1][4] = -det4_0123_0234 * invDet;

	mat[2][0] =  det4_1234_0134 * invDet;
	mat[2][1] = -det4_0234_0134 * invDet;
	mat[2][2] =  det4_0134_0134 * invDet;
	mat[2][3] = -det4_0124_0134 * invDet;
	mat[2][4] =  det4_0123_0134 * invDet;

	mat[3][0] = -det4_1234_0124 * invDet;
	mat[3][1] =  det4_0234_0124 * invDet;
	mat[3][2] = -det4_0134_0124 * invDet;
	mat[3][3] =  det4_0124_0124 * invDet;
	mat[3][4] = -det4_0123_0124 * invDet;

	mat[4][0] =  det4_1234_0123 * invDet;
	mat[4][1] = -det4_0234_0123 * invDet;
	mat[4][2] =  det4_0134_0123 * invDet;
	mat[4][3] = -det4_0124_0123 * invDet;
	mat[4][4] =  det4_0123_0123 * invDet;

	return true;
#elif 0
	// 5*28 = 140 multiplications
	//			5 divisions
	float *mat = reinterpret_cast<float *>(this);
	float s;
	double d, di;

	di = mat[0];
	s = di;
	mat[0] = d = 1.0f / di;
	mat[1] *= d;
	mat[2] *= d;
	mat[3] *= d;
	mat[4] *= d;
	d = -d;
	mat[5] *= d;
	mat[10] *= d;
	mat[15] *= d;
	mat[20] *= d;
	d = mat[5] * di;
	mat[6] += mat[1] * d;
	mat[7] += mat[2] * d;
	mat[8] += mat[3] * d;
	mat[9] += mat[4] * d;
	d = mat[10] * di;
	mat[11] += mat[1] * d;
	mat[12] += mat[2] * d;
	mat[13] += mat[3] * d;
	mat[14] += mat[4] * d;
	d = mat[15] * di;
	mat[16] += mat[1] * d;
	mat[17] += mat[2] * d;
	mat[18] += mat[3] * d;
	mat[19] += mat[4] * d;
	d = mat[20] * di;
	mat[21] += mat[1] * d;
	mat[22] += mat[2] * d;
	mat[23] += mat[3] * d;
	mat[24] += mat[4] * d;
	di = mat[6];
	s *= di;
	mat[6] = d = 1.0f / di;
	mat[5] *= d;
	mat[7] *= d;
	mat[8] *= d;
	mat[9] *= d;
	d = -d;
	mat[1] *= d;
	mat[11] *= d;
	mat[16] *= d;
	mat[21] *= d;
	d = mat[1] * di;
	mat[0] += mat[5] * d;
	mat[2] += mat[7] * d;
	mat[3] += mat[8] * d;
	mat[4] += mat[9] * d;
	d = mat[11] * di;
	mat[10] += mat[5] * d;
	mat[12] += mat[7] * d;
	mat[13] += mat[8] * d;
	mat[14] += mat[9] * d;
	d = mat[16] * di;
	mat[15] += mat[5] * d;
	mat[17] += mat[7] * d;
	mat[18] += mat[8] * d;
	mat[19] += mat[9] * d;
	d = mat[21] * di;
	mat[20] += mat[5] * d;
	mat[22] += mat[7] * d;
	mat[23] += mat[8] * d;
	mat[24] += mat[9] * d;
	di = mat[12];
	s *= di;
	mat[12] = d = 1.0f / di;
	mat[10] *= d;
	mat[11] *= d;
	mat[13] *= d;
	mat[14] *= d;
	d = -d;
	mat[2] *= d;
	mat[7] *= d;
	mat[17] *= d;
	mat[22] *= d;
	d = mat[2] * di;
	mat[0] += mat[10] * d;
	mat[1] += mat[11] * d;
	mat[3] += mat[13] * d;
	mat[4] += mat[14] * d;
	d = mat[7] * di;
	mat[5] += mat[10] * d;
	mat[6] += mat[11] * d;
	mat[8] += mat[13] * d;
	mat[9] += mat[14] * d;
	d = mat[17] * di;
	mat[15] += mat[10] * d;
	mat[16] += mat[11] * d;
	mat[18] += mat[13] * d;
	mat[19] += mat[14] * d;
	d = mat[22] * di;
	mat[20] += mat[10] * d;
	mat[21] += mat[11] * d;
	mat[23] += mat[13] * d;
	mat[24] += mat[14] * d;
	di = mat[18];
	s *= di;
	mat[18] = d = 1.0f / di;
	mat[15] *= d;
	mat[16] *= d;
	mat[17] *= d;
	mat[19] *= d;
	d = -d;
	mat[3] *= d;
	mat[8] *= d;
	mat[13] *= d;
	mat[23] *= d;
	d = mat[3] * di;
	mat[0] += mat[15] * d;
	mat[1] += mat[16] * d;
	mat[2] += mat[17] * d;
	mat[4] += mat[19] * d;
	d = mat[8] * di;
	mat[5] += mat[15] * d;
	mat[6] += mat[16] * d;
	mat[7] += mat[17] * d;
	mat[9] += mat[19] * d;
	d = mat[13] * di;
	mat[10] += mat[15] * d;
	mat[11] += mat[16] * d;
	mat[12] += mat[17] * d;
	mat[14] += mat[19] * d;
	d = mat[23] * di;
	mat[20] += mat[15] * d;
	mat[21] += mat[16] * d;
	mat[22] += mat[17] * d;
	mat[24] += mat[19] * d;
	di = mat[24];
	s *= di;
	mat[24] = d = 1.0f / di;
	mat[20] *= d;
	mat[21] *= d;
	mat[22] *= d;
	mat[23] *= d;
	d = -d;
	mat[4] *= d;
	mat[9] *= d;
	mat[14] *= d;
	mat[19] *= d;
	d = mat[4] * di;
	mat[0] += mat[20] * d;
	mat[1] += mat[21] * d;
	mat[2] += mat[22] * d;
	mat[3] += mat[23] * d;
	d = mat[9] * di;
	mat[5] += mat[20] * d;
	mat[6] += mat[21] * d;
	mat[7] += mat[22] * d;
	mat[8] += mat[23] * d;
	d = mat[14] * di;
	mat[10] += mat[20] * d;
	mat[11] += mat[21] * d;
	mat[12] += mat[22] * d;
	mat[13] += mat[23] * d;
	d = mat[19] * di;
	mat[15] += mat[20] * d;
	mat[16] += mat[21] * d;
	mat[17] += mat[22] * d;
	mat[18] += mat[23] * d;

	return ( s != 0.0f && !FLOAT_IS_NAN( s ) );
#else
	// 86+30+6 = 122 multiplications
	//	  2*1  =   2 divisions
	CMat3D r0, r1, r2, r3;
	float c0, c1, c2, det, invDet;
	float *mat = reinterpret_cast<float *>(this);

	// r0 = m0.Inverse();	// 3x3
	c0 = mat[1*5+1] * mat[2*5+2] - mat[1*5+2] * mat[2*5+1];
	c1 = mat[1*5+2] * mat[2*5+0] - mat[1*5+0] * mat[2*5+2];
	c2 = mat[1*5+0] * mat[2*5+1] - mat[1*5+1] * mat[2*5+0];

	det = mat[0*5+0] * c0 + mat[0*5+1] * c1 + mat[0*5+2] * c2;

	if ( CMath::fabs( det ) < MATRIX_INVERSE_EPSILON ) {
		return false;
	}

	invDet = 1.0f / det;

	r0[0][0] = c0 * invDet;
	r0[0][1] = ( mat[0*5+2] * mat[2*5+1] - mat[0*5+1] * mat[2*5+2] ) * invDet;
	r0[0][2] = ( mat[0*5+1] * mat[1*5+2] - mat[0*5+2] * mat[1*5+1] ) * invDet;
	r0[1][0] = c1 * invDet;
	r0[1][1] = ( mat[0*5+0] * mat[2*5+2] - mat[0*5+2] * mat[2*5+0] ) * invDet;
	r0[1][2] = ( mat[0*5+2] * mat[1*5+0] - mat[0*5+0] * mat[1*5+2] ) * invDet;
	r0[2][0] = c2 * invDet;
	r0[2][1] = ( mat[0*5+1] * mat[2*5+0] - mat[0*5+0] * mat[2*5+1] ) * invDet;
	r0[2][2] = ( mat[0*5+0] * mat[1*5+1] - mat[0*5+1] * mat[1*5+0] ) * invDet;

	// r1 = r0 * m1;		// 3x2 = 3x3 * 3x2
	r1[0][0] = r0[0][0] * mat[0*5+3] + r0[0][1] * mat[1*5+3] + r0[0][2] * mat[2*5+3];
	r1[0][1] = r0[0][0] * mat[0*5+4] + r0[0][1] * mat[1*5+4] + r0[0][2] * mat[2*5+4];
	r1[1][0] = r0[1][0] * mat[0*5+3] + r0[1][1] * mat[1*5+3] + r0[1][2] * mat[2*5+3];
	r1[1][1] = r0[1][0] * mat[0*5+4] + r0[1][1] * mat[1*5+4] + r0[1][2] * mat[2*5+4];
	r1[2][0] = r0[2][0] * mat[0*5+3] + r0[2][1] * mat[1*5+3] + r0[2][2] * mat[2*5+3];
	r1[2][1] = r0[2][0] * mat[0*5+4] + r0[2][1] * mat[1*5+4] + r0[2][2] * mat[2*5+4];

	// r2 = m2 * r1;		// 2x2 = 2x3 * 3x2
	r2[0][0] = mat[3*5+0] * r1[0][0] + mat[3*5+1] * r1[1][0] + mat[3*5+2] * r1[2][0];
	r2[0][1] = mat[3*5+0] * r1[0][1] + mat[3*5+1] * r1[1][1] + mat[3*5+2] * r1[2][1];
	r2[1][0] = mat[4*5+0] * r1[0][0] + mat[4*5+1] * r1[1][0] + mat[4*5+2] * r1[2][0];
	r2[1][1] = mat[4*5+0] * r1[0][1] + mat[4*5+1] * r1[1][1] + mat[4*5+2] * r1[2][1];

	// r3 = r2 - m3;		// 2x2 = 2x2 - 2x2
	r3[0][0] = r2[0][0] - mat[3*5+3];
	r3[0][1] = r2[0][1] - mat[3*5+4];
	r3[1][0] = r2[1][0] - mat[4*5+3];
	r3[1][1] = r2[1][1] - mat[4*5+4];

	// r3.InverseSelf();	// 2x2
	det = r3[0][0] * r3[1][1] - r3[0][1] * r3[1][0];

	if ( CMath::fabs( det ) < MATRIX_INVERSE_EPSILON ) {
		return false;
	}

	invDet = 1.0f / det;

	c0 = r3[0][0];
	r3[0][0] =   r3[1][1] * invDet;
	r3[0][1] = - r3[0][1] * invDet;
	r3[1][0] = - r3[1][0] * invDet;
	r3[1][1] =   c0 * invDet;

	// r2 = m2 * r0;		// 2x3 = 2x3 * 3x3
	r2[0][0] = mat[3*5+0] * r0[0][0] + mat[3*5+1] * r0[1][0] + mat[3*5+2] * r0[2][0];
	r2[0][1] = mat[3*5+0] * r0[0][1] + mat[3*5+1] * r0[1][1] + mat[3*5+2] * r0[2][1];
	r2[0][2] = mat[3*5+0] * r0[0][2] + mat[3*5+1] * r0[1][2] + mat[3*5+2] * r0[2][2];
	r2[1][0] = mat[4*5+0] * r0[0][0] + mat[4*5+1] * r0[1][0] + mat[4*5+2] * r0[2][0];
	r2[1][1] = mat[4*5+0] * r0[0][1] + mat[4*5+1] * r0[1][1] + mat[4*5+2] * r0[2][1];
	r2[1][2] = mat[4*5+0] * r0[0][2] + mat[4*5+1] * r0[1][2] + mat[4*5+2] * r0[2][2];

	// m2 = r3 * r2;		// 2x3 = 2x2 * 2x3
	mat[3*5+0] = r3[0][0] * r2[0][0] + r3[0][1] * r2[1][0];
	mat[3*5+1] = r3[0][0] * r2[0][1] + r3[0][1] * r2[1][1];
	mat[3*5+2] = r3[0][0] * r2[0][2] + r3[0][1] * r2[1][2];
	mat[4*5+0] = r3[1][0] * r2[0][0] + r3[1][1] * r2[1][0];
	mat[4*5+1] = r3[1][0] * r2[0][1] + r3[1][1] * r2[1][1];
	mat[4*5+2] = r3[1][0] * r2[0][2] + r3[1][1] * r2[1][2];

	// m0 = r0 - r1 * m2;	// 3x3 = 3x3 - 3x2 * 2x3
	mat[0*5+0] = r0[0][0] - r1[0][0] * mat[3*5+0] - r1[0][1] * mat[4*5+0];
	mat[0*5+1] = r0[0][1] - r1[0][0] * mat[3*5+1] - r1[0][1] * mat[4*5+1];
	mat[0*5+2] = r0[0][2] - r1[0][0] * mat[3*5+2] - r1[0][1] * mat[4*5+2];
	mat[1*5+0] = r0[1][0] - r1[1][0] * mat[3*5+0] - r1[1][1] * mat[4*5+0];
	mat[1*5+1] = r0[1][1] - r1[1][0] * mat[3*5+1] - r1[1][1] * mat[4*5+1];
	mat[1*5+2] = r0[1][2] - r1[1][0] * mat[3*5+2] - r1[1][1] * mat[4*5+2];
	mat[2*5+0] = r0[2][0] - r1[2][0] * mat[3*5+0] - r1[2][1] * mat[4*5+0];
	mat[2*5+1] = r0[2][1] - r1[2][0] * mat[3*5+1] - r1[2][1] * mat[4*5+1];
	mat[2*5+2] = r0[2][2] - r1[2][0] * mat[3*5+2] - r1[2][1] * mat[4*5+2];

	// m1 = r1 * r3;		// 3x2 = 3x2 * 2x2
	mat[0*5+3] = r1[0][0] * r3[0][0] + r1[0][1] * r3[1][0];
	mat[0*5+4] = r1[0][0] * r3[0][1] + r1[0][1] * r3[1][1];
	mat[1*5+3] = r1[1][0] * r3[0][0] + r1[1][1] * r3[1][0];
	mat[1*5+4] = r1[1][0] * r3[0][1] + r1[1][1] * r3[1][1];
	mat[2*5+3] = r1[2][0] * r3[0][0] + r1[2][1] * r3[1][0];
	mat[2*5+4] = r1[2][0] * r3[0][1] + r1[2][1] * r3[1][1];

	// m3 = -r3;			// 2x2 = - 2x2
	mat[3*5+3] = -r3[0][0];
	mat[3*5+4] = -r3[0][1];
	mat[4*5+3] = -r3[1][0];
	mat[4*5+4] = -r3[1][1];

	return true;
#endif
}

/*
=============
CMat5D::toString
=============
*/
const char *CMat5D::toString( int precision ) const {
	return CMyString::floatArraytoString( toFloatPtr(), getDimension(), precision );
}


//===============================================================
//
//	CMat6D
//
//===============================================================

CMat6D mat6_zero( CVec6D( 0, 0, 0, 0, 0, 0 ), CVec6D( 0, 0, 0, 0, 0, 0 ), CVec6D( 0, 0, 0, 0, 0, 0 ), CVec6D( 0, 0, 0, 0, 0, 0 ), CVec6D( 0, 0, 0, 0, 0, 0 ), CVec6D( 0, 0, 0, 0, 0, 0 ) );
CMat6D mat6_identity( CVec6D( 1, 0, 0, 0, 0, 0 ), CVec6D( 0, 1, 0, 0, 0, 0 ), CVec6D( 0, 0, 1, 0, 0, 0 ), CVec6D( 0, 0, 0, 1, 0, 0 ), CVec6D( 0, 0, 0, 0, 1, 0 ), CVec6D( 0, 0, 0, 0, 0, 1 ) );

/*
============
CMat6D::transpose
============
*/
CMat6D CMat6D::transpose() const {
	CMat6D	transpose;
	int		i, j;

	for( i = 0; i < 6; i++ ) {
		for( j = 0; j < 6; j++ ) {
			transpose[ i ][ j ] = mat[ j ][ i ];
        }
	}
	return transpose;
}

/*
============
CMat6D::transposeSelf
============
*/
CMat6D &CMat6D::transposeSelf() {
	float	temp;
	int		i, j;

	for( i = 0; i < 6; i++ ) {
		for( j = i + 1; j < 6; j++ ) {
			temp = mat[ i ][ j ];
			mat[ i ][ j ] = mat[ j ][ i ];
			mat[ j ][ i ] = temp;
        }
	}
	return *this;
}

/*
============
CMat6D::Determinant
============
*/
float CMat6D::determinant() const {

	// 2x2 sub-determinants required to calculate 6x6 determinant
	float det2_45_01 = mat[4][0] * mat[5][1] - mat[4][1] * mat[5][0];
	float det2_45_02 = mat[4][0] * mat[5][2] - mat[4][2] * mat[5][0];
	float det2_45_03 = mat[4][0] * mat[5][3] - mat[4][3] * mat[5][0];
	float det2_45_04 = mat[4][0] * mat[5][4] - mat[4][4] * mat[5][0];
	float det2_45_05 = mat[4][0] * mat[5][5] - mat[4][5] * mat[5][0];
	float det2_45_12 = mat[4][1] * mat[5][2] - mat[4][2] * mat[5][1];
	float det2_45_13 = mat[4][1] * mat[5][3] - mat[4][3] * mat[5][1];
	float det2_45_14 = mat[4][1] * mat[5][4] - mat[4][4] * mat[5][1];
	float det2_45_15 = mat[4][1] * mat[5][5] - mat[4][5] * mat[5][1];
	float det2_45_23 = mat[4][2] * mat[5][3] - mat[4][3] * mat[5][2];
	float det2_45_24 = mat[4][2] * mat[5][4] - mat[4][4] * mat[5][2];
	float det2_45_25 = mat[4][2] * mat[5][5] - mat[4][5] * mat[5][2];
	float det2_45_34 = mat[4][3] * mat[5][4] - mat[4][4] * mat[5][3];
	float det2_45_35 = mat[4][3] * mat[5][5] - mat[4][5] * mat[5][3];
	float det2_45_45 = mat[4][4] * mat[5][5] - mat[4][5] * mat[5][4];

	// 3x3 sub-determinants required to calculate 6x6 determinant
	float det3_345_012 = mat[3][0] * det2_45_12 - mat[3][1] * det2_45_02 + mat[3][2] * det2_45_01;
	float det3_345_013 = mat[3][0] * det2_45_13 - mat[3][1] * det2_45_03 + mat[3][3] * det2_45_01;
	float det3_345_014 = mat[3][0] * det2_45_14 - mat[3][1] * det2_45_04 + mat[3][4] * det2_45_01;
	float det3_345_015 = mat[3][0] * det2_45_15 - mat[3][1] * det2_45_05 + mat[3][5] * det2_45_01;
	float det3_345_023 = mat[3][0] * det2_45_23 - mat[3][2] * det2_45_03 + mat[3][3] * det2_45_02;
	float det3_345_024 = mat[3][0] * det2_45_24 - mat[3][2] * det2_45_04 + mat[3][4] * det2_45_02;
	float det3_345_025 = mat[3][0] * det2_45_25 - mat[3][2] * det2_45_05 + mat[3][5] * det2_45_02;
	float det3_345_034 = mat[3][0] * det2_45_34 - mat[3][3] * det2_45_04 + mat[3][4] * det2_45_03;
	float det3_345_035 = mat[3][0] * det2_45_35 - mat[3][3] * det2_45_05 + mat[3][5] * det2_45_03;
	float det3_345_045 = mat[3][0] * det2_45_45 - mat[3][4] * det2_45_05 + mat[3][5] * det2_45_04;
	float det3_345_123 = mat[3][1] * det2_45_23 - mat[3][2] * det2_45_13 + mat[3][3] * det2_45_12;
	float det3_345_124 = mat[3][1] * det2_45_24 - mat[3][2] * det2_45_14 + mat[3][4] * det2_45_12;
	float det3_345_125 = mat[3][1] * det2_45_25 - mat[3][2] * det2_45_15 + mat[3][5] * det2_45_12;
	float det3_345_134 = mat[3][1] * det2_45_34 - mat[3][3] * det2_45_14 + mat[3][4] * det2_45_13;
	float det3_345_135 = mat[3][1] * det2_45_35 - mat[3][3] * det2_45_15 + mat[3][5] * det2_45_13;
	float det3_345_145 = mat[3][1] * det2_45_45 - mat[3][4] * det2_45_15 + mat[3][5] * det2_45_14;
	float det3_345_234 = mat[3][2] * det2_45_34 - mat[3][3] * det2_45_24 + mat[3][4] * det2_45_23;
	float det3_345_235 = mat[3][2] * det2_45_35 - mat[3][3] * det2_45_25 + mat[3][5] * det2_45_23;
	float det3_345_245 = mat[3][2] * det2_45_45 - mat[3][4] * det2_45_25 + mat[3][5] * det2_45_24;
	float det3_345_345 = mat[3][3] * det2_45_45 - mat[3][4] * det2_45_35 + mat[3][5] * det2_45_34;

	// 4x4 sub-determinants required to calculate 6x6 determinant
	float det4_2345_0123 = mat[2][0] * det3_345_123 - mat[2][1] * det3_345_023 + mat[2][2] * det3_345_013 - mat[2][3] * det3_345_012;
	float det4_2345_0124 = mat[2][0] * det3_345_124 - mat[2][1] * det3_345_024 + mat[2][2] * det3_345_014 - mat[2][4] * det3_345_012;
	float det4_2345_0125 = mat[2][0] * det3_345_125 - mat[2][1] * det3_345_025 + mat[2][2] * det3_345_015 - mat[2][5] * det3_345_012;
	float det4_2345_0134 = mat[2][0] * det3_345_134 - mat[2][1] * det3_345_034 + mat[2][3] * det3_345_014 - mat[2][4] * det3_345_013;
	float det4_2345_0135 = mat[2][0] * det3_345_135 - mat[2][1] * det3_345_035 + mat[2][3] * det3_345_015 - mat[2][5] * det3_345_013;
	float det4_2345_0145 = mat[2][0] * det3_345_145 - mat[2][1] * det3_345_045 + mat[2][4] * det3_345_015 - mat[2][5] * det3_345_014;
	float det4_2345_0234 = mat[2][0] * det3_345_234 - mat[2][2] * det3_345_034 + mat[2][3] * det3_345_024 - mat[2][4] * det3_345_023;
	float det4_2345_0235 = mat[2][0] * det3_345_235 - mat[2][2] * det3_345_035 + mat[2][3] * det3_345_025 - mat[2][5] * det3_345_023;
	float det4_2345_0245 = mat[2][0] * det3_345_245 - mat[2][2] * det3_345_045 + mat[2][4] * det3_345_025 - mat[2][5] * det3_345_024;
	float det4_2345_0345 = mat[2][0] * det3_345_345 - mat[2][3] * det3_345_045 + mat[2][4] * det3_345_035 - mat[2][5] * det3_345_034;
	float det4_2345_1234 = mat[2][1] * det3_345_234 - mat[2][2] * det3_345_134 + mat[2][3] * det3_345_124 - mat[2][4] * det3_345_123;
	float det4_2345_1235 = mat[2][1] * det3_345_235 - mat[2][2] * det3_345_135 + mat[2][3] * det3_345_125 - mat[2][5] * det3_345_123;
	float det4_2345_1245 = mat[2][1] * det3_345_245 - mat[2][2] * det3_345_145 + mat[2][4] * det3_345_125 - mat[2][5] * det3_345_124;
	float det4_2345_1345 = mat[2][1] * det3_345_345 - mat[2][3] * det3_345_145 + mat[2][4] * det3_345_135 - mat[2][5] * det3_345_134;
	float det4_2345_2345 = mat[2][2] * det3_345_345 - mat[2][3] * det3_345_245 + mat[2][4] * det3_345_235 - mat[2][5] * det3_345_234;

	// 5x5 sub-determinants required to calculate 6x6 determinant
	float det5_12345_01234 = mat[1][0] * det4_2345_1234 - mat[1][1] * det4_2345_0234 + mat[1][2] * det4_2345_0134 - mat[1][3] * det4_2345_0124 + mat[1][4] * det4_2345_0123;
	float det5_12345_01235 = mat[1][0] * det4_2345_1235 - mat[1][1] * det4_2345_0235 + mat[1][2] * det4_2345_0135 - mat[1][3] * det4_2345_0125 + mat[1][5] * det4_2345_0123;
	float det5_12345_01245 = mat[1][0] * det4_2345_1245 - mat[1][1] * det4_2345_0245 + mat[1][2] * det4_2345_0145 - mat[1][4] * det4_2345_0125 + mat[1][5] * det4_2345_0124;
	float det5_12345_01345 = mat[1][0] * det4_2345_1345 - mat[1][1] * det4_2345_0345 + mat[1][3] * det4_2345_0145 - mat[1][4] * det4_2345_0135 + mat[1][5] * det4_2345_0134;
	float det5_12345_02345 = mat[1][0] * det4_2345_2345 - mat[1][2] * det4_2345_0345 + mat[1][3] * det4_2345_0245 - mat[1][4] * det4_2345_0235 + mat[1][5] * det4_2345_0234;
	float det5_12345_12345 = mat[1][1] * det4_2345_2345 - mat[1][2] * det4_2345_1345 + mat[1][3] * det4_2345_1245 - mat[1][4] * det4_2345_1235 + mat[1][5] * det4_2345_1234;

	// determinant of 6x6 matrix
	return	mat[0][0] * det5_12345_12345 - mat[0][1] * det5_12345_02345 + mat[0][2] * det5_12345_01345 -
			mat[0][3] * det5_12345_01245 + mat[0][4] * det5_12345_01235 - mat[0][5] * det5_12345_01234;
}

/*
============
CMat6D::InverseSelf
============
*/
bool CMat6D::inverseSelf() {
	// 810+6+36 = 852 multiplications
	//				1 division
	double det, invDet;

	// 2x2 sub-determinants required to calculate 6x6 determinant
	float det2_45_01 = mat[4][0] * mat[5][1] - mat[4][1] * mat[5][0];
	float det2_45_02 = mat[4][0] * mat[5][2] - mat[4][2] * mat[5][0];
	float det2_45_03 = mat[4][0] * mat[5][3] - mat[4][3] * mat[5][0];
	float det2_45_04 = mat[4][0] * mat[5][4] - mat[4][4] * mat[5][0];
	float det2_45_05 = mat[4][0] * mat[5][5] - mat[4][5] * mat[5][0];
	float det2_45_12 = mat[4][1] * mat[5][2] - mat[4][2] * mat[5][1];
	float det2_45_13 = mat[4][1] * mat[5][3] - mat[4][3] * mat[5][1];
	float det2_45_14 = mat[4][1] * mat[5][4] - mat[4][4] * mat[5][1];
	float det2_45_15 = mat[4][1] * mat[5][5] - mat[4][5] * mat[5][1];
	float det2_45_23 = mat[4][2] * mat[5][3] - mat[4][3] * mat[5][2];
	float det2_45_24 = mat[4][2] * mat[5][4] - mat[4][4] * mat[5][2];
	float det2_45_25 = mat[4][2] * mat[5][5] - mat[4][5] * mat[5][2];
	float det2_45_34 = mat[4][3] * mat[5][4] - mat[4][4] * mat[5][3];
	float det2_45_35 = mat[4][3] * mat[5][5] - mat[4][5] * mat[5][3];
	float det2_45_45 = mat[4][4] * mat[5][5] - mat[4][5] * mat[5][4];

	// 3x3 sub-determinants required to calculate 6x6 determinant
	float det3_345_012 = mat[3][0] * det2_45_12 - mat[3][1] * det2_45_02 + mat[3][2] * det2_45_01;
	float det3_345_013 = mat[3][0] * det2_45_13 - mat[3][1] * det2_45_03 + mat[3][3] * det2_45_01;
	float det3_345_014 = mat[3][0] * det2_45_14 - mat[3][1] * det2_45_04 + mat[3][4] * det2_45_01;
	float det3_345_015 = mat[3][0] * det2_45_15 - mat[3][1] * det2_45_05 + mat[3][5] * det2_45_01;
	float det3_345_023 = mat[3][0] * det2_45_23 - mat[3][2] * det2_45_03 + mat[3][3] * det2_45_02;
	float det3_345_024 = mat[3][0] * det2_45_24 - mat[3][2] * det2_45_04 + mat[3][4] * det2_45_02;
	float det3_345_025 = mat[3][0] * det2_45_25 - mat[3][2] * det2_45_05 + mat[3][5] * det2_45_02;
	float det3_345_034 = mat[3][0] * det2_45_34 - mat[3][3] * det2_45_04 + mat[3][4] * det2_45_03;
	float det3_345_035 = mat[3][0] * det2_45_35 - mat[3][3] * det2_45_05 + mat[3][5] * det2_45_03;
	float det3_345_045 = mat[3][0] * det2_45_45 - mat[3][4] * det2_45_05 + mat[3][5] * det2_45_04;
	float det3_345_123 = mat[3][1] * det2_45_23 - mat[3][2] * det2_45_13 + mat[3][3] * det2_45_12;
	float det3_345_124 = mat[3][1] * det2_45_24 - mat[3][2] * det2_45_14 + mat[3][4] * det2_45_12;
	float det3_345_125 = mat[3][1] * det2_45_25 - mat[3][2] * det2_45_15 + mat[3][5] * det2_45_12;
	float det3_345_134 = mat[3][1] * det2_45_34 - mat[3][3] * det2_45_14 + mat[3][4] * det2_45_13;
	float det3_345_135 = mat[3][1] * det2_45_35 - mat[3][3] * det2_45_15 + mat[3][5] * det2_45_13;
	float det3_345_145 = mat[3][1] * det2_45_45 - mat[3][4] * det2_45_15 + mat[3][5] * det2_45_14;
	float det3_345_234 = mat[3][2] * det2_45_34 - mat[3][3] * det2_45_24 + mat[3][4] * det2_45_23;
	float det3_345_235 = mat[3][2] * det2_45_35 - mat[3][3] * det2_45_25 + mat[3][5] * det2_45_23;
	float det3_345_245 = mat[3][2] * det2_45_45 - mat[3][4] * det2_45_25 + mat[3][5] * det2_45_24;
	float det3_345_345 = mat[3][3] * det2_45_45 - mat[3][4] * det2_45_35 + mat[3][5] * det2_45_34;

	// 4x4 sub-determinants required to calculate 6x6 determinant
	float det4_2345_0123 = mat[2][0] * det3_345_123 - mat[2][1] * det3_345_023 + mat[2][2] * det3_345_013 - mat[2][3] * det3_345_012;
	float det4_2345_0124 = mat[2][0] * det3_345_124 - mat[2][1] * det3_345_024 + mat[2][2] * det3_345_014 - mat[2][4] * det3_345_012;
	float det4_2345_0125 = mat[2][0] * det3_345_125 - mat[2][1] * det3_345_025 + mat[2][2] * det3_345_015 - mat[2][5] * det3_345_012;
	float det4_2345_0134 = mat[2][0] * det3_345_134 - mat[2][1] * det3_345_034 + mat[2][3] * det3_345_014 - mat[2][4] * det3_345_013;
	float det4_2345_0135 = mat[2][0] * det3_345_135 - mat[2][1] * det3_345_035 + mat[2][3] * det3_345_015 - mat[2][5] * det3_345_013;
	float det4_2345_0145 = mat[2][0] * det3_345_145 - mat[2][1] * det3_345_045 + mat[2][4] * det3_345_015 - mat[2][5] * det3_345_014;
	float det4_2345_0234 = mat[2][0] * det3_345_234 - mat[2][2] * det3_345_034 + mat[2][3] * det3_345_024 - mat[2][4] * det3_345_023;
	float det4_2345_0235 = mat[2][0] * det3_345_235 - mat[2][2] * det3_345_035 + mat[2][3] * det3_345_025 - mat[2][5] * det3_345_023;
	float det4_2345_0245 = mat[2][0] * det3_345_245 - mat[2][2] * det3_345_045 + mat[2][4] * det3_345_025 - mat[2][5] * det3_345_024;
	float det4_2345_0345 = mat[2][0] * det3_345_345 - mat[2][3] * det3_345_045 + mat[2][4] * det3_345_035 - mat[2][5] * det3_345_034;
	float det4_2345_1234 = mat[2][1] * det3_345_234 - mat[2][2] * det3_345_134 + mat[2][3] * det3_345_124 - mat[2][4] * det3_345_123;
	float det4_2345_1235 = mat[2][1] * det3_345_235 - mat[2][2] * det3_345_135 + mat[2][3] * det3_345_125 - mat[2][5] * det3_345_123;
	float det4_2345_1245 = mat[2][1] * det3_345_245 - mat[2][2] * det3_345_145 + mat[2][4] * det3_345_125 - mat[2][5] * det3_345_124;
	float det4_2345_1345 = mat[2][1] * det3_345_345 - mat[2][3] * det3_345_145 + mat[2][4] * det3_345_135 - mat[2][5] * det3_345_134;
	float det4_2345_2345 = mat[2][2] * det3_345_345 - mat[2][3] * det3_345_245 + mat[2][4] * det3_345_235 - mat[2][5] * det3_345_234;

	// 5x5 sub-determinants required to calculate 6x6 determinant
	float det5_12345_01234 = mat[1][0] * det4_2345_1234 - mat[1][1] * det4_2345_0234 + mat[1][2] * det4_2345_0134 - mat[1][3] * det4_2345_0124 + mat[1][4] * det4_2345_0123;
	float det5_12345_01235 = mat[1][0] * det4_2345_1235 - mat[1][1] * det4_2345_0235 + mat[1][2] * det4_2345_0135 - mat[1][3] * det4_2345_0125 + mat[1][5] * det4_2345_0123;
	float det5_12345_01245 = mat[1][0] * det4_2345_1245 - mat[1][1] * det4_2345_0245 + mat[1][2] * det4_2345_0145 - mat[1][4] * det4_2345_0125 + mat[1][5] * det4_2345_0124;
	float det5_12345_01345 = mat[1][0] * det4_2345_1345 - mat[1][1] * det4_2345_0345 + mat[1][3] * det4_2345_0145 - mat[1][4] * det4_2345_0135 + mat[1][5] * det4_2345_0134;
	float det5_12345_02345 = mat[1][0] * det4_2345_2345 - mat[1][2] * det4_2345_0345 + mat[1][3] * det4_2345_0245 - mat[1][4] * det4_2345_0235 + mat[1][5] * det4_2345_0234;
	float det5_12345_12345 = mat[1][1] * det4_2345_2345 - mat[1][2] * det4_2345_1345 + mat[1][3] * det4_2345_1245 - mat[1][4] * det4_2345_1235 + mat[1][5] * det4_2345_1234;

	// determinant of 6x6 matrix
	det = mat[0][0] * det5_12345_12345 - mat[0][1] * det5_12345_02345 + mat[0][2] * det5_12345_01345 -
				mat[0][3] * det5_12345_01245 + mat[0][4] * det5_12345_01235 - mat[0][5] * det5_12345_01234;

	if ( CMath::fabs( det ) < MATRIX_INVERSE_EPSILON ) {
		return false;
	}

	invDet = 1.0f / det;

	// remaining 2x2 sub-determinants
	float det2_34_01 = mat[3][0] * mat[4][1] - mat[3][1] * mat[4][0];
	float det2_34_02 = mat[3][0] * mat[4][2] - mat[3][2] * mat[4][0];
	float det2_34_03 = mat[3][0] * mat[4][3] - mat[3][3] * mat[4][0];
	float det2_34_04 = mat[3][0] * mat[4][4] - mat[3][4] * mat[4][0];
	float det2_34_05 = mat[3][0] * mat[4][5] - mat[3][5] * mat[4][0];
	float det2_34_12 = mat[3][1] * mat[4][2] - mat[3][2] * mat[4][1];
	float det2_34_13 = mat[3][1] * mat[4][3] - mat[3][3] * mat[4][1];
	float det2_34_14 = mat[3][1] * mat[4][4] - mat[3][4] * mat[4][1];
	float det2_34_15 = mat[3][1] * mat[4][5] - mat[3][5] * mat[4][1];
	float det2_34_23 = mat[3][2] * mat[4][3] - mat[3][3] * mat[4][2];
	float det2_34_24 = mat[3][2] * mat[4][4] - mat[3][4] * mat[4][2];
	float det2_34_25 = mat[3][2] * mat[4][5] - mat[3][5] * mat[4][2];
	float det2_34_34 = mat[3][3] * mat[4][4] - mat[3][4] * mat[4][3];
	float det2_34_35 = mat[3][3] * mat[4][5] - mat[3][5] * mat[4][3];
	float det2_34_45 = mat[3][4] * mat[4][5] - mat[3][5] * mat[4][4];
	float det2_35_01 = mat[3][0] * mat[5][1] - mat[3][1] * mat[5][0];
	float det2_35_02 = mat[3][0] * mat[5][2] - mat[3][2] * mat[5][0];
	float det2_35_03 = mat[3][0] * mat[5][3] - mat[3][3] * mat[5][0];
	float det2_35_04 = mat[3][0] * mat[5][4] - mat[3][4] * mat[5][0];
	float det2_35_05 = mat[3][0] * mat[5][5] - mat[3][5] * mat[5][0];
	float det2_35_12 = mat[3][1] * mat[5][2] - mat[3][2] * mat[5][1];
	float det2_35_13 = mat[3][1] * mat[5][3] - mat[3][3] * mat[5][1];
	float det2_35_14 = mat[3][1] * mat[5][4] - mat[3][4] * mat[5][1];
	float det2_35_15 = mat[3][1] * mat[5][5] - mat[3][5] * mat[5][1];
	float det2_35_23 = mat[3][2] * mat[5][3] - mat[3][3] * mat[5][2];
	float det2_35_24 = mat[3][2] * mat[5][4] - mat[3][4] * mat[5][2];
	float det2_35_25 = mat[3][2] * mat[5][5] - mat[3][5] * mat[5][2];
	float det2_35_34 = mat[3][3] * mat[5][4] - mat[3][4] * mat[5][3];
	float det2_35_35 = mat[3][3] * mat[5][5] - mat[3][5] * mat[5][3];
	float det2_35_45 = mat[3][4] * mat[5][5] - mat[3][5] * mat[5][4];

	// remaining 3x3 sub-determinants
	float det3_234_012 = mat[2][0] * det2_34_12 - mat[2][1] * det2_34_02 + mat[2][2] * det2_34_01;
	float det3_234_013 = mat[2][0] * det2_34_13 - mat[2][1] * det2_34_03 + mat[2][3] * det2_34_01;
	float det3_234_014 = mat[2][0] * det2_34_14 - mat[2][1] * det2_34_04 + mat[2][4] * det2_34_01;
	float det3_234_015 = mat[2][0] * det2_34_15 - mat[2][1] * det2_34_05 + mat[2][5] * det2_34_01;
	float det3_234_023 = mat[2][0] * det2_34_23 - mat[2][2] * det2_34_03 + mat[2][3] * det2_34_02;
	float det3_234_024 = mat[2][0] * det2_34_24 - mat[2][2] * det2_34_04 + mat[2][4] * det2_34_02;
	float det3_234_025 = mat[2][0] * det2_34_25 - mat[2][2] * det2_34_05 + mat[2][5] * det2_34_02;
	float det3_234_034 = mat[2][0] * det2_34_34 - mat[2][3] * det2_34_04 + mat[2][4] * det2_34_03;
	float det3_234_035 = mat[2][0] * det2_34_35 - mat[2][3] * det2_34_05 + mat[2][5] * det2_34_03;
	float det3_234_045 = mat[2][0] * det2_34_45 - mat[2][4] * det2_34_05 + mat[2][5] * det2_34_04;
	float det3_234_123 = mat[2][1] * det2_34_23 - mat[2][2] * det2_34_13 + mat[2][3] * det2_34_12;
	float det3_234_124 = mat[2][1] * det2_34_24 - mat[2][2] * det2_34_14 + mat[2][4] * det2_34_12;
	float det3_234_125 = mat[2][1] * det2_34_25 - mat[2][2] * det2_34_15 + mat[2][5] * det2_34_12;
	float det3_234_134 = mat[2][1] * det2_34_34 - mat[2][3] * det2_34_14 + mat[2][4] * det2_34_13;
	float det3_234_135 = mat[2][1] * det2_34_35 - mat[2][3] * det2_34_15 + mat[2][5] * det2_34_13;
	float det3_234_145 = mat[2][1] * det2_34_45 - mat[2][4] * det2_34_15 + mat[2][5] * det2_34_14;
	float det3_234_234 = mat[2][2] * det2_34_34 - mat[2][3] * det2_34_24 + mat[2][4] * det2_34_23;
	float det3_234_235 = mat[2][2] * det2_34_35 - mat[2][3] * det2_34_25 + mat[2][5] * det2_34_23;
	float det3_234_245 = mat[2][2] * det2_34_45 - mat[2][4] * det2_34_25 + mat[2][5] * det2_34_24;
	float det3_234_345 = mat[2][3] * det2_34_45 - mat[2][4] * det2_34_35 + mat[2][5] * det2_34_34;
	float det3_235_012 = mat[2][0] * det2_35_12 - mat[2][1] * det2_35_02 + mat[2][2] * det2_35_01;
	float det3_235_013 = mat[2][0] * det2_35_13 - mat[2][1] * det2_35_03 + mat[2][3] * det2_35_01;
	float det3_235_014 = mat[2][0] * det2_35_14 - mat[2][1] * det2_35_04 + mat[2][4] * det2_35_01;
	float det3_235_015 = mat[2][0] * det2_35_15 - mat[2][1] * det2_35_05 + mat[2][5] * det2_35_01;
	float det3_235_023 = mat[2][0] * det2_35_23 - mat[2][2] * det2_35_03 + mat[2][3] * det2_35_02;
	float det3_235_024 = mat[2][0] * det2_35_24 - mat[2][2] * det2_35_04 + mat[2][4] * det2_35_02;
	float det3_235_025 = mat[2][0] * det2_35_25 - mat[2][2] * det2_35_05 + mat[2][5] * det2_35_02;
	float det3_235_034 = mat[2][0] * det2_35_34 - mat[2][3] * det2_35_04 + mat[2][4] * det2_35_03;
	float det3_235_035 = mat[2][0] * det2_35_35 - mat[2][3] * det2_35_05 + mat[2][5] * det2_35_03;
	float det3_235_045 = mat[2][0] * det2_35_45 - mat[2][4] * det2_35_05 + mat[2][5] * det2_35_04;
	float det3_235_123 = mat[2][1] * det2_35_23 - mat[2][2] * det2_35_13 + mat[2][3] * det2_35_12;
	float det3_235_124 = mat[2][1] * det2_35_24 - mat[2][2] * det2_35_14 + mat[2][4] * det2_35_12;
	float det3_235_125 = mat[2][1] * det2_35_25 - mat[2][2] * det2_35_15 + mat[2][5] * det2_35_12;
	float det3_235_134 = mat[2][1] * det2_35_34 - mat[2][3] * det2_35_14 + mat[2][4] * det2_35_13;
	float det3_235_135 = mat[2][1] * det2_35_35 - mat[2][3] * det2_35_15 + mat[2][5] * det2_35_13;
	float det3_235_145 = mat[2][1] * det2_35_45 - mat[2][4] * det2_35_15 + mat[2][5] * det2_35_14;
	float det3_235_234 = mat[2][2] * det2_35_34 - mat[2][3] * det2_35_24 + mat[2][4] * det2_35_23;
	float det3_235_235 = mat[2][2] * det2_35_35 - mat[2][3] * det2_35_25 + mat[2][5] * det2_35_23;
	float det3_235_245 = mat[2][2] * det2_35_45 - mat[2][4] * det2_35_25 + mat[2][5] * det2_35_24;
	float det3_235_345 = mat[2][3] * det2_35_45 - mat[2][4] * det2_35_35 + mat[2][5] * det2_35_34;
	float det3_245_012 = mat[2][0] * det2_45_12 - mat[2][1] * det2_45_02 + mat[2][2] * det2_45_01;
	float det3_245_013 = mat[2][0] * det2_45_13 - mat[2][1] * det2_45_03 + mat[2][3] * det2_45_01;
	float det3_245_014 = mat[2][0] * det2_45_14 - mat[2][1] * det2_45_04 + mat[2][4] * det2_45_01;
	float det3_245_015 = mat[2][0] * det2_45_15 - mat[2][1] * det2_45_05 + mat[2][5] * det2_45_01;
	float det3_245_023 = mat[2][0] * det2_45_23 - mat[2][2] * det2_45_03 + mat[2][3] * det2_45_02;
	float det3_245_024 = mat[2][0] * det2_45_24 - mat[2][2] * det2_45_04 + mat[2][4] * det2_45_02;
	float det3_245_025 = mat[2][0] * det2_45_25 - mat[2][2] * det2_45_05 + mat[2][5] * det2_45_02;
	float det3_245_034 = mat[2][0] * det2_45_34 - mat[2][3] * det2_45_04 + mat[2][4] * det2_45_03;
	float det3_245_035 = mat[2][0] * det2_45_35 - mat[2][3] * det2_45_05 + mat[2][5] * det2_45_03;
	float det3_245_045 = mat[2][0] * det2_45_45 - mat[2][4] * det2_45_05 + mat[2][5] * det2_45_04;
	float det3_245_123 = mat[2][1] * det2_45_23 - mat[2][2] * det2_45_13 + mat[2][3] * det2_45_12;
	float det3_245_124 = mat[2][1] * det2_45_24 - mat[2][2] * det2_45_14 + mat[2][4] * det2_45_12;
	float det3_245_125 = mat[2][1] * det2_45_25 - mat[2][2] * det2_45_15 + mat[2][5] * det2_45_12;
	float det3_245_134 = mat[2][1] * det2_45_34 - mat[2][3] * det2_45_14 + mat[2][4] * det2_45_13;
	float det3_245_135 = mat[2][1] * det2_45_35 - mat[2][3] * det2_45_15 + mat[2][5] * det2_45_13;
	float det3_245_145 = mat[2][1] * det2_45_45 - mat[2][4] * det2_45_15 + mat[2][5] * det2_45_14;
	float det3_245_234 = mat[2][2] * det2_45_34 - mat[2][3] * det2_45_24 + mat[2][4] * det2_45_23;
	float det3_245_235 = mat[2][2] * det2_45_35 - mat[2][3] * det2_45_25 + mat[2][5] * det2_45_23;
	float det3_245_245 = mat[2][2] * det2_45_45 - mat[2][4] * det2_45_25 + mat[2][5] * det2_45_24;
	float det3_245_345 = mat[2][3] * det2_45_45 - mat[2][4] * det2_45_35 + mat[2][5] * det2_45_34;

	// remaining 4x4 sub-determinants
	float det4_1234_0123 = mat[1][0] * det3_234_123 - mat[1][1] * det3_234_023 + mat[1][2] * det3_234_013 - mat[1][3] * det3_234_012;
	float det4_1234_0124 = mat[1][0] * det3_234_124 - mat[1][1] * det3_234_024 + mat[1][2] * det3_234_014 - mat[1][4] * det3_234_012;
	float det4_1234_0125 = mat[1][0] * det3_234_125 - mat[1][1] * det3_234_025 + mat[1][2] * det3_234_015 - mat[1][5] * det3_234_012;
	float det4_1234_0134 = mat[1][0] * det3_234_134 - mat[1][1] * det3_234_034 + mat[1][3] * det3_234_014 - mat[1][4] * det3_234_013;
	float det4_1234_0135 = mat[1][0] * det3_234_135 - mat[1][1] * det3_234_035 + mat[1][3] * det3_234_015 - mat[1][5] * det3_234_013;
	float det4_1234_0145 = mat[1][0] * det3_234_145 - mat[1][1] * det3_234_045 + mat[1][4] * det3_234_015 - mat[1][5] * det3_234_014;
	float det4_1234_0234 = mat[1][0] * det3_234_234 - mat[1][2] * det3_234_034 + mat[1][3] * det3_234_024 - mat[1][4] * det3_234_023;
	float det4_1234_0235 = mat[1][0] * det3_234_235 - mat[1][2] * det3_234_035 + mat[1][3] * det3_234_025 - mat[1][5] * det3_234_023;
	float det4_1234_0245 = mat[1][0] * det3_234_245 - mat[1][2] * det3_234_045 + mat[1][4] * det3_234_025 - mat[1][5] * det3_234_024;
	float det4_1234_0345 = mat[1][0] * det3_234_345 - mat[1][3] * det3_234_045 + mat[1][4] * det3_234_035 - mat[1][5] * det3_234_034;
	float det4_1234_1234 = mat[1][1] * det3_234_234 - mat[1][2] * det3_234_134 + mat[1][3] * det3_234_124 - mat[1][4] * det3_234_123;
	float det4_1234_1235 = mat[1][1] * det3_234_235 - mat[1][2] * det3_234_135 + mat[1][3] * det3_234_125 - mat[1][5] * det3_234_123;
	float det4_1234_1245 = mat[1][1] * det3_234_245 - mat[1][2] * det3_234_145 + mat[1][4] * det3_234_125 - mat[1][5] * det3_234_124;
	float det4_1234_1345 = mat[1][1] * det3_234_345 - mat[1][3] * det3_234_145 + mat[1][4] * det3_234_135 - mat[1][5] * det3_234_134;
	float det4_1234_2345 = mat[1][2] * det3_234_345 - mat[1][3] * det3_234_245 + mat[1][4] * det3_234_235 - mat[1][5] * det3_234_234;
	float det4_1235_0123 = mat[1][0] * det3_235_123 - mat[1][1] * det3_235_023 + mat[1][2] * det3_235_013 - mat[1][3] * det3_235_012;
	float det4_1235_0124 = mat[1][0] * det3_235_124 - mat[1][1] * det3_235_024 + mat[1][2] * det3_235_014 - mat[1][4] * det3_235_012;
	float det4_1235_0125 = mat[1][0] * det3_235_125 - mat[1][1] * det3_235_025 + mat[1][2] * det3_235_015 - mat[1][5] * det3_235_012;
	float det4_1235_0134 = mat[1][0] * det3_235_134 - mat[1][1] * det3_235_034 + mat[1][3] * det3_235_014 - mat[1][4] * det3_235_013;
	float det4_1235_0135 = mat[1][0] * det3_235_135 - mat[1][1] * det3_235_035 + mat[1][3] * det3_235_015 - mat[1][5] * det3_235_013;
	float det4_1235_0145 = mat[1][0] * det3_235_145 - mat[1][1] * det3_235_045 + mat[1][4] * det3_235_015 - mat[1][5] * det3_235_014;
	float det4_1235_0234 = mat[1][0] * det3_235_234 - mat[1][2] * det3_235_034 + mat[1][3] * det3_235_024 - mat[1][4] * det3_235_023;
	float det4_1235_0235 = mat[1][0] * det3_235_235 - mat[1][2] * det3_235_035 + mat[1][3] * det3_235_025 - mat[1][5] * det3_235_023;
	float det4_1235_0245 = mat[1][0] * det3_235_245 - mat[1][2] * det3_235_045 + mat[1][4] * det3_235_025 - mat[1][5] * det3_235_024;
	float det4_1235_0345 = mat[1][0] * det3_235_345 - mat[1][3] * det3_235_045 + mat[1][4] * det3_235_035 - mat[1][5] * det3_235_034;
	float det4_1235_1234 = mat[1][1] * det3_235_234 - mat[1][2] * det3_235_134 + mat[1][3] * det3_235_124 - mat[1][4] * det3_235_123;
	float det4_1235_1235 = mat[1][1] * det3_235_235 - mat[1][2] * det3_235_135 + mat[1][3] * det3_235_125 - mat[1][5] * det3_235_123;
	float det4_1235_1245 = mat[1][1] * det3_235_245 - mat[1][2] * det3_235_145 + mat[1][4] * det3_235_125 - mat[1][5] * det3_235_124;
	float det4_1235_1345 = mat[1][1] * det3_235_345 - mat[1][3] * det3_235_145 + mat[1][4] * det3_235_135 - mat[1][5] * det3_235_134;
	float det4_1235_2345 = mat[1][2] * det3_235_345 - mat[1][3] * det3_235_245 + mat[1][4] * det3_235_235 - mat[1][5] * det3_235_234;
	float det4_1245_0123 = mat[1][0] * det3_245_123 - mat[1][1] * det3_245_023 + mat[1][2] * det3_245_013 - mat[1][3] * det3_245_012;
	float det4_1245_0124 = mat[1][0] * det3_245_124 - mat[1][1] * det3_245_024 + mat[1][2] * det3_245_014 - mat[1][4] * det3_245_012;
	float det4_1245_0125 = mat[1][0] * det3_245_125 - mat[1][1] * det3_245_025 + mat[1][2] * det3_245_015 - mat[1][5] * det3_245_012;
	float det4_1245_0134 = mat[1][0] * det3_245_134 - mat[1][1] * det3_245_034 + mat[1][3] * det3_245_014 - mat[1][4] * det3_245_013;
	float det4_1245_0135 = mat[1][0] * det3_245_135 - mat[1][1] * det3_245_035 + mat[1][3] * det3_245_015 - mat[1][5] * det3_245_013;
	float det4_1245_0145 = mat[1][0] * det3_245_145 - mat[1][1] * det3_245_045 + mat[1][4] * det3_245_015 - mat[1][5] * det3_245_014;
	float det4_1245_0234 = mat[1][0] * det3_245_234 - mat[1][2] * det3_245_034 + mat[1][3] * det3_245_024 - mat[1][4] * det3_245_023;
	float det4_1245_0235 = mat[1][0] * det3_245_235 - mat[1][2] * det3_245_035 + mat[1][3] * det3_245_025 - mat[1][5] * det3_245_023;
	float det4_1245_0245 = mat[1][0] * det3_245_245 - mat[1][2] * det3_245_045 + mat[1][4] * det3_245_025 - mat[1][5] * det3_245_024;
	float det4_1245_0345 = mat[1][0] * det3_245_345 - mat[1][3] * det3_245_045 + mat[1][4] * det3_245_035 - mat[1][5] * det3_245_034;
	float det4_1245_1234 = mat[1][1] * det3_245_234 - mat[1][2] * det3_245_134 + mat[1][3] * det3_245_124 - mat[1][4] * det3_245_123;
	float det4_1245_1235 = mat[1][1] * det3_245_235 - mat[1][2] * det3_245_135 + mat[1][3] * det3_245_125 - mat[1][5] * det3_245_123;
	float det4_1245_1245 = mat[1][1] * det3_245_245 - mat[1][2] * det3_245_145 + mat[1][4] * det3_245_125 - mat[1][5] * det3_245_124;
	float det4_1245_1345 = mat[1][1] * det3_245_345 - mat[1][3] * det3_245_145 + mat[1][4] * det3_245_135 - mat[1][5] * det3_245_134;
	float det4_1245_2345 = mat[1][2] * det3_245_345 - mat[1][3] * det3_245_245 + mat[1][4] * det3_245_235 - mat[1][5] * det3_245_234;
	float det4_1345_0123 = mat[1][0] * det3_345_123 - mat[1][1] * det3_345_023 + mat[1][2] * det3_345_013 - mat[1][3] * det3_345_012;
	float det4_1345_0124 = mat[1][0] * det3_345_124 - mat[1][1] * det3_345_024 + mat[1][2] * det3_345_014 - mat[1][4] * det3_345_012;
	float det4_1345_0125 = mat[1][0] * det3_345_125 - mat[1][1] * det3_345_025 + mat[1][2] * det3_345_015 - mat[1][5] * det3_345_012;
	float det4_1345_0134 = mat[1][0] * det3_345_134 - mat[1][1] * det3_345_034 + mat[1][3] * det3_345_014 - mat[1][4] * det3_345_013;
	float det4_1345_0135 = mat[1][0] * det3_345_135 - mat[1][1] * det3_345_035 + mat[1][3] * det3_345_015 - mat[1][5] * det3_345_013;
	float det4_1345_0145 = mat[1][0] * det3_345_145 - mat[1][1] * det3_345_045 + mat[1][4] * det3_345_015 - mat[1][5] * det3_345_014;
	float det4_1345_0234 = mat[1][0] * det3_345_234 - mat[1][2] * det3_345_034 + mat[1][3] * det3_345_024 - mat[1][4] * det3_345_023;
	float det4_1345_0235 = mat[1][0] * det3_345_235 - mat[1][2] * det3_345_035 + mat[1][3] * det3_345_025 - mat[1][5] * det3_345_023;
	float det4_1345_0245 = mat[1][0] * det3_345_245 - mat[1][2] * det3_345_045 + mat[1][4] * det3_345_025 - mat[1][5] * det3_345_024;
	float det4_1345_0345 = mat[1][0] * det3_345_345 - mat[1][3] * det3_345_045 + mat[1][4] * det3_345_035 - mat[1][5] * det3_345_034;
	float det4_1345_1234 = mat[1][1] * det3_345_234 - mat[1][2] * det3_345_134 + mat[1][3] * det3_345_124 - mat[1][4] * det3_345_123;
	float det4_1345_1235 = mat[1][1] * det3_345_235 - mat[1][2] * det3_345_135 + mat[1][3] * det3_345_125 - mat[1][5] * det3_345_123;
	float det4_1345_1245 = mat[1][1] * det3_345_245 - mat[1][2] * det3_345_145 + mat[1][4] * det3_345_125 - mat[1][5] * det3_345_124;
	float det4_1345_1345 = mat[1][1] * det3_345_345 - mat[1][3] * det3_345_145 + mat[1][4] * det3_345_135 - mat[1][5] * det3_345_134;
	float det4_1345_2345 = mat[1][2] * det3_345_345 - mat[1][3] * det3_345_245 + mat[1][4] * det3_345_235 - mat[1][5] * det3_345_234;

	// remaining 5x5 sub-determinants
	float det5_01234_01234 = mat[0][0] * det4_1234_1234 - mat[0][1] * det4_1234_0234 + mat[0][2] * det4_1234_0134 - mat[0][3] * det4_1234_0124 + mat[0][4] * det4_1234_0123;
	float det5_01234_01235 = mat[0][0] * det4_1234_1235 - mat[0][1] * det4_1234_0235 + mat[0][2] * det4_1234_0135 - mat[0][3] * det4_1234_0125 + mat[0][5] * det4_1234_0123;
	float det5_01234_01245 = mat[0][0] * det4_1234_1245 - mat[0][1] * det4_1234_0245 + mat[0][2] * det4_1234_0145 - mat[0][4] * det4_1234_0125 + mat[0][5] * det4_1234_0124;
	float det5_01234_01345 = mat[0][0] * det4_1234_1345 - mat[0][1] * det4_1234_0345 + mat[0][3] * det4_1234_0145 - mat[0][4] * det4_1234_0135 + mat[0][5] * det4_1234_0134;
	float det5_01234_02345 = mat[0][0] * det4_1234_2345 - mat[0][2] * det4_1234_0345 + mat[0][3] * det4_1234_0245 - mat[0][4] * det4_1234_0235 + mat[0][5] * det4_1234_0234;
	float det5_01234_12345 = mat[0][1] * det4_1234_2345 - mat[0][2] * det4_1234_1345 + mat[0][3] * det4_1234_1245 - mat[0][4] * det4_1234_1235 + mat[0][5] * det4_1234_1234;
	float det5_01235_01234 = mat[0][0] * det4_1235_1234 - mat[0][1] * det4_1235_0234 + mat[0][2] * det4_1235_0134 - mat[0][3] * det4_1235_0124 + mat[0][4] * det4_1235_0123;
	float det5_01235_01235 = mat[0][0] * det4_1235_1235 - mat[0][1] * det4_1235_0235 + mat[0][2] * det4_1235_0135 - mat[0][3] * det4_1235_0125 + mat[0][5] * det4_1235_0123;
	float det5_01235_01245 = mat[0][0] * det4_1235_1245 - mat[0][1] * det4_1235_0245 + mat[0][2] * det4_1235_0145 - mat[0][4] * det4_1235_0125 + mat[0][5] * det4_1235_0124;
	float det5_01235_01345 = mat[0][0] * det4_1235_1345 - mat[0][1] * det4_1235_0345 + mat[0][3] * det4_1235_0145 - mat[0][4] * det4_1235_0135 + mat[0][5] * det4_1235_0134;
	float det5_01235_02345 = mat[0][0] * det4_1235_2345 - mat[0][2] * det4_1235_0345 + mat[0][3] * det4_1235_0245 - mat[0][4] * det4_1235_0235 + mat[0][5] * det4_1235_0234;
	float det5_01235_12345 = mat[0][1] * det4_1235_2345 - mat[0][2] * det4_1235_1345 + mat[0][3] * det4_1235_1245 - mat[0][4] * det4_1235_1235 + mat[0][5] * det4_1235_1234;
	float det5_01245_01234 = mat[0][0] * det4_1245_1234 - mat[0][1] * det4_1245_0234 + mat[0][2] * det4_1245_0134 - mat[0][3] * det4_1245_0124 + mat[0][4] * det4_1245_0123;
	float det5_01245_01235 = mat[0][0] * det4_1245_1235 - mat[0][1] * det4_1245_0235 + mat[0][2] * det4_1245_0135 - mat[0][3] * det4_1245_0125 + mat[0][5] * det4_1245_0123;
	float det5_01245_01245 = mat[0][0] * det4_1245_1245 - mat[0][1] * det4_1245_0245 + mat[0][2] * det4_1245_0145 - mat[0][4] * det4_1245_0125 + mat[0][5] * det4_1245_0124;
	float det5_01245_01345 = mat[0][0] * det4_1245_1345 - mat[0][1] * det4_1245_0345 + mat[0][3] * det4_1245_0145 - mat[0][4] * det4_1245_0135 + mat[0][5] * det4_1245_0134;
	float det5_01245_02345 = mat[0][0] * det4_1245_2345 - mat[0][2] * det4_1245_0345 + mat[0][3] * det4_1245_0245 - mat[0][4] * det4_1245_0235 + mat[0][5] * det4_1245_0234;
	float det5_01245_12345 = mat[0][1] * det4_1245_2345 - mat[0][2] * det4_1245_1345 + mat[0][3] * det4_1245_1245 - mat[0][4] * det4_1245_1235 + mat[0][5] * det4_1245_1234;
	float det5_01345_01234 = mat[0][0] * det4_1345_1234 - mat[0][1] * det4_1345_0234 + mat[0][2] * det4_1345_0134 - mat[0][3] * det4_1345_0124 + mat[0][4] * det4_1345_0123;
	float det5_01345_01235 = mat[0][0] * det4_1345_1235 - mat[0][1] * det4_1345_0235 + mat[0][2] * det4_1345_0135 - mat[0][3] * det4_1345_0125 + mat[0][5] * det4_1345_0123;
	float det5_01345_01245 = mat[0][0] * det4_1345_1245 - mat[0][1] * det4_1345_0245 + mat[0][2] * det4_1345_0145 - mat[0][4] * det4_1345_0125 + mat[0][5] * det4_1345_0124;
	float det5_01345_01345 = mat[0][0] * det4_1345_1345 - mat[0][1] * det4_1345_0345 + mat[0][3] * det4_1345_0145 - mat[0][4] * det4_1345_0135 + mat[0][5] * det4_1345_0134;
	float det5_01345_02345 = mat[0][0] * det4_1345_2345 - mat[0][2] * det4_1345_0345 + mat[0][3] * det4_1345_0245 - mat[0][4] * det4_1345_0235 + mat[0][5] * det4_1345_0234;
	float det5_01345_12345 = mat[0][1] * det4_1345_2345 - mat[0][2] * det4_1345_1345 + mat[0][3] * det4_1345_1245 - mat[0][4] * det4_1345_1235 + mat[0][5] * det4_1345_1234;
	float det5_02345_01234 = mat[0][0] * det4_2345_1234 - mat[0][1] * det4_2345_0234 + mat[0][2] * det4_2345_0134 - mat[0][3] * det4_2345_0124 + mat[0][4] * det4_2345_0123;
	float det5_02345_01235 = mat[0][0] * det4_2345_1235 - mat[0][1] * det4_2345_0235 + mat[0][2] * det4_2345_0135 - mat[0][3] * det4_2345_0125 + mat[0][5] * det4_2345_0123;
	float det5_02345_01245 = mat[0][0] * det4_2345_1245 - mat[0][1] * det4_2345_0245 + mat[0][2] * det4_2345_0145 - mat[0][4] * det4_2345_0125 + mat[0][5] * det4_2345_0124;
	float det5_02345_01345 = mat[0][0] * det4_2345_1345 - mat[0][1] * det4_2345_0345 + mat[0][3] * det4_2345_0145 - mat[0][4] * det4_2345_0135 + mat[0][5] * det4_2345_0134;
	float det5_02345_02345 = mat[0][0] * det4_2345_2345 - mat[0][2] * det4_2345_0345 + mat[0][3] * det4_2345_0245 - mat[0][4] * det4_2345_0235 + mat[0][5] * det4_2345_0234;
	float det5_02345_12345 = mat[0][1] * det4_2345_2345 - mat[0][2] * det4_2345_1345 + mat[0][3] * det4_2345_1245 - mat[0][4] * det4_2345_1235 + mat[0][5] * det4_2345_1234;

	mat[0][0] =  det5_12345_12345 * invDet;
	mat[0][1] = -det5_02345_12345 * invDet;
	mat[0][2] =  det5_01345_12345 * invDet;
	mat[0][3] = -det5_01245_12345 * invDet;
	mat[0][4] =  det5_01235_12345 * invDet;
	mat[0][5] = -det5_01234_12345 * invDet;

	mat[1][0] = -det5_12345_02345 * invDet;
	mat[1][1] =  det5_02345_02345 * invDet;
	mat[1][2] = -det5_01345_02345 * invDet;
	mat[1][3] =  det5_01245_02345 * invDet;
	mat[1][4] = -det5_01235_02345 * invDet;
	mat[1][5] =  det5_01234_02345 * invDet;

	mat[2][0] =  det5_12345_01345 * invDet;
	mat[2][1] = -det5_02345_01345 * invDet;
	mat[2][2] =  det5_01345_01345 * invDet;
	mat[2][3] = -det5_01245_01345 * invDet;
	mat[2][4] =  det5_01235_01345 * invDet;
	mat[2][5] = -det5_01234_01345 * invDet;

	mat[3][0] = -det5_12345_01245 * invDet;
	mat[3][1] =  det5_02345_01245 * invDet;
	mat[3][2] = -det5_01345_01245 * invDet;
	mat[3][3] =  det5_01245_01245 * invDet;
	mat[3][4] = -det5_01235_01245 * invDet;
	mat[3][5] =  det5_01234_01245 * invDet;

	mat[4][0] =  det5_12345_01235 * invDet;
	mat[4][1] = -det5_02345_01235 * invDet;
	mat[4][2] =  det5_01345_01235 * invDet;
	mat[4][3] = -det5_01245_01235 * invDet;
	mat[4][4] =  det5_01235_01235 * invDet;
	mat[4][5] = -det5_01234_01235 * invDet;

	mat[5][0] = -det5_12345_01234 * invDet;
	mat[5][1] =  det5_02345_01234 * invDet;
	mat[5][2] = -det5_01345_01234 * invDet;
	mat[5][3] =  det5_01245_01234 * invDet;
	mat[5][4] = -det5_01235_01234 * invDet;
	mat[5][5] =  det5_01234_01234 * invDet;

	return true;
}

/*
============
CMat6D::InverseFastSelf
============
*/
bool CMat6D::inverseFastSelf() {
#if 0
	// 810+6+36 = 852 multiplications
	//				1 division
	double det, invDet;

	// 2x2 sub-determinants required to calculate 6x6 determinant
	float det2_45_01 = mat[4][0] * mat[5][1] - mat[4][1] * mat[5][0];
	float det2_45_02 = mat[4][0] * mat[5][2] - mat[4][2] * mat[5][0];
	float det2_45_03 = mat[4][0] * mat[5][3] - mat[4][3] * mat[5][0];
	float det2_45_04 = mat[4][0] * mat[5][4] - mat[4][4] * mat[5][0];
	float det2_45_05 = mat[4][0] * mat[5][5] - mat[4][5] * mat[5][0];
	float det2_45_12 = mat[4][1] * mat[5][2] - mat[4][2] * mat[5][1];
	float det2_45_13 = mat[4][1] * mat[5][3] - mat[4][3] * mat[5][1];
	float det2_45_14 = mat[4][1] * mat[5][4] - mat[4][4] * mat[5][1];
	float det2_45_15 = mat[4][1] * mat[5][5] - mat[4][5] * mat[5][1];
	float det2_45_23 = mat[4][2] * mat[5][3] - mat[4][3] * mat[5][2];
	float det2_45_24 = mat[4][2] * mat[5][4] - mat[4][4] * mat[5][2];
	float det2_45_25 = mat[4][2] * mat[5][5] - mat[4][5] * mat[5][2];
	float det2_45_34 = mat[4][3] * mat[5][4] - mat[4][4] * mat[5][3];
	float det2_45_35 = mat[4][3] * mat[5][5] - mat[4][5] * mat[5][3];
	float det2_45_45 = mat[4][4] * mat[5][5] - mat[4][5] * mat[5][4];

	// 3x3 sub-determinants required to calculate 6x6 determinant
	float det3_345_012 = mat[3][0] * det2_45_12 - mat[3][1] * det2_45_02 + mat[3][2] * det2_45_01;
	float det3_345_013 = mat[3][0] * det2_45_13 - mat[3][1] * det2_45_03 + mat[3][3] * det2_45_01;
	float det3_345_014 = mat[3][0] * det2_45_14 - mat[3][1] * det2_45_04 + mat[3][4] * det2_45_01;
	float det3_345_015 = mat[3][0] * det2_45_15 - mat[3][1] * det2_45_05 + mat[3][5] * det2_45_01;
	float det3_345_023 = mat[3][0] * det2_45_23 - mat[3][2] * det2_45_03 + mat[3][3] * det2_45_02;
	float det3_345_024 = mat[3][0] * det2_45_24 - mat[3][2] * det2_45_04 + mat[3][4] * det2_45_02;
	float det3_345_025 = mat[3][0] * det2_45_25 - mat[3][2] * det2_45_05 + mat[3][5] * det2_45_02;
	float det3_345_034 = mat[3][0] * det2_45_34 - mat[3][3] * det2_45_04 + mat[3][4] * det2_45_03;
	float det3_345_035 = mat[3][0] * det2_45_35 - mat[3][3] * det2_45_05 + mat[3][5] * det2_45_03;
	float det3_345_045 = mat[3][0] * det2_45_45 - mat[3][4] * det2_45_05 + mat[3][5] * det2_45_04;
	float det3_345_123 = mat[3][1] * det2_45_23 - mat[3][2] * det2_45_13 + mat[3][3] * det2_45_12;
	float det3_345_124 = mat[3][1] * det2_45_24 - mat[3][2] * det2_45_14 + mat[3][4] * det2_45_12;
	float det3_345_125 = mat[3][1] * det2_45_25 - mat[3][2] * det2_45_15 + mat[3][5] * det2_45_12;
	float det3_345_134 = mat[3][1] * det2_45_34 - mat[3][3] * det2_45_14 + mat[3][4] * det2_45_13;
	float det3_345_135 = mat[3][1] * det2_45_35 - mat[3][3] * det2_45_15 + mat[3][5] * det2_45_13;
	float det3_345_145 = mat[3][1] * det2_45_45 - mat[3][4] * det2_45_15 + mat[3][5] * det2_45_14;
	float det3_345_234 = mat[3][2] * det2_45_34 - mat[3][3] * det2_45_24 + mat[3][4] * det2_45_23;
	float det3_345_235 = mat[3][2] * det2_45_35 - mat[3][3] * det2_45_25 + mat[3][5] * det2_45_23;
	float det3_345_245 = mat[3][2] * det2_45_45 - mat[3][4] * det2_45_25 + mat[3][5] * det2_45_24;
	float det3_345_345 = mat[3][3] * det2_45_45 - mat[3][4] * det2_45_35 + mat[3][5] * det2_45_34;

	// 4x4 sub-determinants required to calculate 6x6 determinant
	float det4_2345_0123 = mat[2][0] * det3_345_123 - mat[2][1] * det3_345_023 + mat[2][2] * det3_345_013 - mat[2][3] * det3_345_012;
	float det4_2345_0124 = mat[2][0] * det3_345_124 - mat[2][1] * det3_345_024 + mat[2][2] * det3_345_014 - mat[2][4] * det3_345_012;
	float det4_2345_0125 = mat[2][0] * det3_345_125 - mat[2][1] * det3_345_025 + mat[2][2] * det3_345_015 - mat[2][5] * det3_345_012;
	float det4_2345_0134 = mat[2][0] * det3_345_134 - mat[2][1] * det3_345_034 + mat[2][3] * det3_345_014 - mat[2][4] * det3_345_013;
	float det4_2345_0135 = mat[2][0] * det3_345_135 - mat[2][1] * det3_345_035 + mat[2][3] * det3_345_015 - mat[2][5] * det3_345_013;
	float det4_2345_0145 = mat[2][0] * det3_345_145 - mat[2][1] * det3_345_045 + mat[2][4] * det3_345_015 - mat[2][5] * det3_345_014;
	float det4_2345_0234 = mat[2][0] * det3_345_234 - mat[2][2] * det3_345_034 + mat[2][3] * det3_345_024 - mat[2][4] * det3_345_023;
	float det4_2345_0235 = mat[2][0] * det3_345_235 - mat[2][2] * det3_345_035 + mat[2][3] * det3_345_025 - mat[2][5] * det3_345_023;
	float det4_2345_0245 = mat[2][0] * det3_345_245 - mat[2][2] * det3_345_045 + mat[2][4] * det3_345_025 - mat[2][5] * det3_345_024;
	float det4_2345_0345 = mat[2][0] * det3_345_345 - mat[2][3] * det3_345_045 + mat[2][4] * det3_345_035 - mat[2][5] * det3_345_034;
	float det4_2345_1234 = mat[2][1] * det3_345_234 - mat[2][2] * det3_345_134 + mat[2][3] * det3_345_124 - mat[2][4] * det3_345_123;
	float det4_2345_1235 = mat[2][1] * det3_345_235 - mat[2][2] * det3_345_135 + mat[2][3] * det3_345_125 - mat[2][5] * det3_345_123;
	float det4_2345_1245 = mat[2][1] * det3_345_245 - mat[2][2] * det3_345_145 + mat[2][4] * det3_345_125 - mat[2][5] * det3_345_124;
	float det4_2345_1345 = mat[2][1] * det3_345_345 - mat[2][3] * det3_345_145 + mat[2][4] * det3_345_135 - mat[2][5] * det3_345_134;
	float det4_2345_2345 = mat[2][2] * det3_345_345 - mat[2][3] * det3_345_245 + mat[2][4] * det3_345_235 - mat[2][5] * det3_345_234;

	// 5x5 sub-determinants required to calculate 6x6 determinant
	float det5_12345_01234 = mat[1][0] * det4_2345_1234 - mat[1][1] * det4_2345_0234 + mat[1][2] * det4_2345_0134 - mat[1][3] * det4_2345_0124 + mat[1][4] * det4_2345_0123;
	float det5_12345_01235 = mat[1][0] * det4_2345_1235 - mat[1][1] * det4_2345_0235 + mat[1][2] * det4_2345_0135 - mat[1][3] * det4_2345_0125 + mat[1][5] * det4_2345_0123;
	float det5_12345_01245 = mat[1][0] * det4_2345_1245 - mat[1][1] * det4_2345_0245 + mat[1][2] * det4_2345_0145 - mat[1][4] * det4_2345_0125 + mat[1][5] * det4_2345_0124;
	float det5_12345_01345 = mat[1][0] * det4_2345_1345 - mat[1][1] * det4_2345_0345 + mat[1][3] * det4_2345_0145 - mat[1][4] * det4_2345_0135 + mat[1][5] * det4_2345_0134;
	float det5_12345_02345 = mat[1][0] * det4_2345_2345 - mat[1][2] * det4_2345_0345 + mat[1][3] * det4_2345_0245 - mat[1][4] * det4_2345_0235 + mat[1][5] * det4_2345_0234;
	float det5_12345_12345 = mat[1][1] * det4_2345_2345 - mat[1][2] * det4_2345_1345 + mat[1][3] * det4_2345_1245 - mat[1][4] * det4_2345_1235 + mat[1][5] * det4_2345_1234;

	// determinant of 6x6 matrix
	det = mat[0][0] * det5_12345_12345 - mat[0][1] * det5_12345_02345 + mat[0][2] * det5_12345_01345 -
				mat[0][3] * det5_12345_01245 + mat[0][4] * det5_12345_01235 - mat[0][5] * det5_12345_01234;

	if ( CMath::fabs( det ) < MATRIX_INVERSE_EPSILON ) {
		return false;
	}

	invDet = 1.0f / det;

	// remaining 2x2 sub-determinants
	float det2_34_01 = mat[3][0] * mat[4][1] - mat[3][1] * mat[4][0];
	float det2_34_02 = mat[3][0] * mat[4][2] - mat[3][2] * mat[4][0];
	float det2_34_03 = mat[3][0] * mat[4][3] - mat[3][3] * mat[4][0];
	float det2_34_04 = mat[3][0] * mat[4][4] - mat[3][4] * mat[4][0];
	float det2_34_05 = mat[3][0] * mat[4][5] - mat[3][5] * mat[4][0];
	float det2_34_12 = mat[3][1] * mat[4][2] - mat[3][2] * mat[4][1];
	float det2_34_13 = mat[3][1] * mat[4][3] - mat[3][3] * mat[4][1];
	float det2_34_14 = mat[3][1] * mat[4][4] - mat[3][4] * mat[4][1];
	float det2_34_15 = mat[3][1] * mat[4][5] - mat[3][5] * mat[4][1];
	float det2_34_23 = mat[3][2] * mat[4][3] - mat[3][3] * mat[4][2];
	float det2_34_24 = mat[3][2] * mat[4][4] - mat[3][4] * mat[4][2];
	float det2_34_25 = mat[3][2] * mat[4][5] - mat[3][5] * mat[4][2];
	float det2_34_34 = mat[3][3] * mat[4][4] - mat[3][4] * mat[4][3];
	float det2_34_35 = mat[3][3] * mat[4][5] - mat[3][5] * mat[4][3];
	float det2_34_45 = mat[3][4] * mat[4][5] - mat[3][5] * mat[4][4];
	float det2_35_01 = mat[3][0] * mat[5][1] - mat[3][1] * mat[5][0];
	float det2_35_02 = mat[3][0] * mat[5][2] - mat[3][2] * mat[5][0];
	float det2_35_03 = mat[3][0] * mat[5][3] - mat[3][3] * mat[5][0];
	float det2_35_04 = mat[3][0] * mat[5][4] - mat[3][4] * mat[5][0];
	float det2_35_05 = mat[3][0] * mat[5][5] - mat[3][5] * mat[5][0];
	float det2_35_12 = mat[3][1] * mat[5][2] - mat[3][2] * mat[5][1];
	float det2_35_13 = mat[3][1] * mat[5][3] - mat[3][3] * mat[5][1];
	float det2_35_14 = mat[3][1] * mat[5][4] - mat[3][4] * mat[5][1];
	float det2_35_15 = mat[3][1] * mat[5][5] - mat[3][5] * mat[5][1];
	float det2_35_23 = mat[3][2] * mat[5][3] - mat[3][3] * mat[5][2];
	float det2_35_24 = mat[3][2] * mat[5][4] - mat[3][4] * mat[5][2];
	float det2_35_25 = mat[3][2] * mat[5][5] - mat[3][5] * mat[5][2];
	float det2_35_34 = mat[3][3] * mat[5][4] - mat[3][4] * mat[5][3];
	float det2_35_35 = mat[3][3] * mat[5][5] - mat[3][5] * mat[5][3];
	float det2_35_45 = mat[3][4] * mat[5][5] - mat[3][5] * mat[5][4];

	// remaining 3x3 sub-determinants
	float det3_234_012 = mat[2][0] * det2_34_12 - mat[2][1] * det2_34_02 + mat[2][2] * det2_34_01;
	float det3_234_013 = mat[2][0] * det2_34_13 - mat[2][1] * det2_34_03 + mat[2][3] * det2_34_01;
	float det3_234_014 = mat[2][0] * det2_34_14 - mat[2][1] * det2_34_04 + mat[2][4] * det2_34_01;
	float det3_234_015 = mat[2][0] * det2_34_15 - mat[2][1] * det2_34_05 + mat[2][5] * det2_34_01;
	float det3_234_023 = mat[2][0] * det2_34_23 - mat[2][2] * det2_34_03 + mat[2][3] * det2_34_02;
	float det3_234_024 = mat[2][0] * det2_34_24 - mat[2][2] * det2_34_04 + mat[2][4] * det2_34_02;
	float det3_234_025 = mat[2][0] * det2_34_25 - mat[2][2] * det2_34_05 + mat[2][5] * det2_34_02;
	float det3_234_034 = mat[2][0] * det2_34_34 - mat[2][3] * det2_34_04 + mat[2][4] * det2_34_03;
	float det3_234_035 = mat[2][0] * det2_34_35 - mat[2][3] * det2_34_05 + mat[2][5] * det2_34_03;
	float det3_234_045 = mat[2][0] * det2_34_45 - mat[2][4] * det2_34_05 + mat[2][5] * det2_34_04;
	float det3_234_123 = mat[2][1] * det2_34_23 - mat[2][2] * det2_34_13 + mat[2][3] * det2_34_12;
	float det3_234_124 = mat[2][1] * det2_34_24 - mat[2][2] * det2_34_14 + mat[2][4] * det2_34_12;
	float det3_234_125 = mat[2][1] * det2_34_25 - mat[2][2] * det2_34_15 + mat[2][5] * det2_34_12;
	float det3_234_134 = mat[2][1] * det2_34_34 - mat[2][3] * det2_34_14 + mat[2][4] * det2_34_13;
	float det3_234_135 = mat[2][1] * det2_34_35 - mat[2][3] * det2_34_15 + mat[2][5] * det2_34_13;
	float det3_234_145 = mat[2][1] * det2_34_45 - mat[2][4] * det2_34_15 + mat[2][5] * det2_34_14;
	float det3_234_234 = mat[2][2] * det2_34_34 - mat[2][3] * det2_34_24 + mat[2][4] * det2_34_23;
	float det3_234_235 = mat[2][2] * det2_34_35 - mat[2][3] * det2_34_25 + mat[2][5] * det2_34_23;
	float det3_234_245 = mat[2][2] * det2_34_45 - mat[2][4] * det2_34_25 + mat[2][5] * det2_34_24;
	float det3_234_345 = mat[2][3] * det2_34_45 - mat[2][4] * det2_34_35 + mat[2][5] * det2_34_34;
	float det3_235_012 = mat[2][0] * det2_35_12 - mat[2][1] * det2_35_02 + mat[2][2] * det2_35_01;
	float det3_235_013 = mat[2][0] * det2_35_13 - mat[2][1] * det2_35_03 + mat[2][3] * det2_35_01;
	float det3_235_014 = mat[2][0] * det2_35_14 - mat[2][1] * det2_35_04 + mat[2][4] * det2_35_01;
	float det3_235_015 = mat[2][0] * det2_35_15 - mat[2][1] * det2_35_05 + mat[2][5] * det2_35_01;
	float det3_235_023 = mat[2][0] * det2_35_23 - mat[2][2] * det2_35_03 + mat[2][3] * det2_35_02;
	float det3_235_024 = mat[2][0] * det2_35_24 - mat[2][2] * det2_35_04 + mat[2][4] * det2_35_02;
	float det3_235_025 = mat[2][0] * det2_35_25 - mat[2][2] * det2_35_05 + mat[2][5] * det2_35_02;
	float det3_235_034 = mat[2][0] * det2_35_34 - mat[2][3] * det2_35_04 + mat[2][4] * det2_35_03;
	float det3_235_035 = mat[2][0] * det2_35_35 - mat[2][3] * det2_35_05 + mat[2][5] * det2_35_03;
	float det3_235_045 = mat[2][0] * det2_35_45 - mat[2][4] * det2_35_05 + mat[2][5] * det2_35_04;
	float det3_235_123 = mat[2][1] * det2_35_23 - mat[2][2] * det2_35_13 + mat[2][3] * det2_35_12;
	float det3_235_124 = mat[2][1] * det2_35_24 - mat[2][2] * det2_35_14 + mat[2][4] * det2_35_12;
	float det3_235_125 = mat[2][1] * det2_35_25 - mat[2][2] * det2_35_15 + mat[2][5] * det2_35_12;
	float det3_235_134 = mat[2][1] * det2_35_34 - mat[2][3] * det2_35_14 + mat[2][4] * det2_35_13;
	float det3_235_135 = mat[2][1] * det2_35_35 - mat[2][3] * det2_35_15 + mat[2][5] * det2_35_13;
	float det3_235_145 = mat[2][1] * det2_35_45 - mat[2][4] * det2_35_15 + mat[2][5] * det2_35_14;
	float det3_235_234 = mat[2][2] * det2_35_34 - mat[2][3] * det2_35_24 + mat[2][4] * det2_35_23;
	float det3_235_235 = mat[2][2] * det2_35_35 - mat[2][3] * det2_35_25 + mat[2][5] * det2_35_23;
	float det3_235_245 = mat[2][2] * det2_35_45 - mat[2][4] * det2_35_25 + mat[2][5] * det2_35_24;
	float det3_235_345 = mat[2][3] * det2_35_45 - mat[2][4] * det2_35_35 + mat[2][5] * det2_35_34;
	float det3_245_012 = mat[2][0] * det2_45_12 - mat[2][1] * det2_45_02 + mat[2][2] * det2_45_01;
	float det3_245_013 = mat[2][0] * det2_45_13 - mat[2][1] * det2_45_03 + mat[2][3] * det2_45_01;
	float det3_245_014 = mat[2][0] * det2_45_14 - mat[2][1] * det2_45_04 + mat[2][4] * det2_45_01;
	float det3_245_015 = mat[2][0] * det2_45_15 - mat[2][1] * det2_45_05 + mat[2][5] * det2_45_01;
	float det3_245_023 = mat[2][0] * det2_45_23 - mat[2][2] * det2_45_03 + mat[2][3] * det2_45_02;
	float det3_245_024 = mat[2][0] * det2_45_24 - mat[2][2] * det2_45_04 + mat[2][4] * det2_45_02;
	float det3_245_025 = mat[2][0] * det2_45_25 - mat[2][2] * det2_45_05 + mat[2][5] * det2_45_02;
	float det3_245_034 = mat[2][0] * det2_45_34 - mat[2][3] * det2_45_04 + mat[2][4] * det2_45_03;
	float det3_245_035 = mat[2][0] * det2_45_35 - mat[2][3] * det2_45_05 + mat[2][5] * det2_45_03;
	float det3_245_045 = mat[2][0] * det2_45_45 - mat[2][4] * det2_45_05 + mat[2][5] * det2_45_04;
	float det3_245_123 = mat[2][1] * det2_45_23 - mat[2][2] * det2_45_13 + mat[2][3] * det2_45_12;
	float det3_245_124 = mat[2][1] * det2_45_24 - mat[2][2] * det2_45_14 + mat[2][4] * det2_45_12;
	float det3_245_125 = mat[2][1] * det2_45_25 - mat[2][2] * det2_45_15 + mat[2][5] * det2_45_12;
	float det3_245_134 = mat[2][1] * det2_45_34 - mat[2][3] * det2_45_14 + mat[2][4] * det2_45_13;
	float det3_245_135 = mat[2][1] * det2_45_35 - mat[2][3] * det2_45_15 + mat[2][5] * det2_45_13;
	float det3_245_145 = mat[2][1] * det2_45_45 - mat[2][4] * det2_45_15 + mat[2][5] * det2_45_14;
	float det3_245_234 = mat[2][2] * det2_45_34 - mat[2][3] * det2_45_24 + mat[2][4] * det2_45_23;
	float det3_245_235 = mat[2][2] * det2_45_35 - mat[2][3] * det2_45_25 + mat[2][5] * det2_45_23;
	float det3_245_245 = mat[2][2] * det2_45_45 - mat[2][4] * det2_45_25 + mat[2][5] * det2_45_24;
	float det3_245_345 = mat[2][3] * det2_45_45 - mat[2][4] * det2_45_35 + mat[2][5] * det2_45_34;

	// remaining 4x4 sub-determinants
	float det4_1234_0123 = mat[1][0] * det3_234_123 - mat[1][1] * det3_234_023 + mat[1][2] * det3_234_013 - mat[1][3] * det3_234_012;
	float det4_1234_0124 = mat[1][0] * det3_234_124 - mat[1][1] * det3_234_024 + mat[1][2] * det3_234_014 - mat[1][4] * det3_234_012;
	float det4_1234_0125 = mat[1][0] * det3_234_125 - mat[1][1] * det3_234_025 + mat[1][2] * det3_234_015 - mat[1][5] * det3_234_012;
	float det4_1234_0134 = mat[1][0] * det3_234_134 - mat[1][1] * det3_234_034 + mat[1][3] * det3_234_014 - mat[1][4] * det3_234_013;
	float det4_1234_0135 = mat[1][0] * det3_234_135 - mat[1][1] * det3_234_035 + mat[1][3] * det3_234_015 - mat[1][5] * det3_234_013;
	float det4_1234_0145 = mat[1][0] * det3_234_145 - mat[1][1] * det3_234_045 + mat[1][4] * det3_234_015 - mat[1][5] * det3_234_014;
	float det4_1234_0234 = mat[1][0] * det3_234_234 - mat[1][2] * det3_234_034 + mat[1][3] * det3_234_024 - mat[1][4] * det3_234_023;
	float det4_1234_0235 = mat[1][0] * det3_234_235 - mat[1][2] * det3_234_035 + mat[1][3] * det3_234_025 - mat[1][5] * det3_234_023;
	float det4_1234_0245 = mat[1][0] * det3_234_245 - mat[1][2] * det3_234_045 + mat[1][4] * det3_234_025 - mat[1][5] * det3_234_024;
	float det4_1234_0345 = mat[1][0] * det3_234_345 - mat[1][3] * det3_234_045 + mat[1][4] * det3_234_035 - mat[1][5] * det3_234_034;
	float det4_1234_1234 = mat[1][1] * det3_234_234 - mat[1][2] * det3_234_134 + mat[1][3] * det3_234_124 - mat[1][4] * det3_234_123;
	float det4_1234_1235 = mat[1][1] * det3_234_235 - mat[1][2] * det3_234_135 + mat[1][3] * det3_234_125 - mat[1][5] * det3_234_123;
	float det4_1234_1245 = mat[1][1] * det3_234_245 - mat[1][2] * det3_234_145 + mat[1][4] * det3_234_125 - mat[1][5] * det3_234_124;
	float det4_1234_1345 = mat[1][1] * det3_234_345 - mat[1][3] * det3_234_145 + mat[1][4] * det3_234_135 - mat[1][5] * det3_234_134;
	float det4_1234_2345 = mat[1][2] * det3_234_345 - mat[1][3] * det3_234_245 + mat[1][4] * det3_234_235 - mat[1][5] * det3_234_234;
	float det4_1235_0123 = mat[1][0] * det3_235_123 - mat[1][1] * det3_235_023 + mat[1][2] * det3_235_013 - mat[1][3] * det3_235_012;
	float det4_1235_0124 = mat[1][0] * det3_235_124 - mat[1][1] * det3_235_024 + mat[1][2] * det3_235_014 - mat[1][4] * det3_235_012;
	float det4_1235_0125 = mat[1][0] * det3_235_125 - mat[1][1] * det3_235_025 + mat[1][2] * det3_235_015 - mat[1][5] * det3_235_012;
	float det4_1235_0134 = mat[1][0] * det3_235_134 - mat[1][1] * det3_235_034 + mat[1][3] * det3_235_014 - mat[1][4] * det3_235_013;
	float det4_1235_0135 = mat[1][0] * det3_235_135 - mat[1][1] * det3_235_035 + mat[1][3] * det3_235_015 - mat[1][5] * det3_235_013;
	float det4_1235_0145 = mat[1][0] * det3_235_145 - mat[1][1] * det3_235_045 + mat[1][4] * det3_235_015 - mat[1][5] * det3_235_014;
	float det4_1235_0234 = mat[1][0] * det3_235_234 - mat[1][2] * det3_235_034 + mat[1][3] * det3_235_024 - mat[1][4] * det3_235_023;
	float det4_1235_0235 = mat[1][0] * det3_235_235 - mat[1][2] * det3_235_035 + mat[1][3] * det3_235_025 - mat[1][5] * det3_235_023;
	float det4_1235_0245 = mat[1][0] * det3_235_245 - mat[1][2] * det3_235_045 + mat[1][4] * det3_235_025 - mat[1][5] * det3_235_024;
	float det4_1235_0345 = mat[1][0] * det3_235_345 - mat[1][3] * det3_235_045 + mat[1][4] * det3_235_035 - mat[1][5] * det3_235_034;
	float det4_1235_1234 = mat[1][1] * det3_235_234 - mat[1][2] * det3_235_134 + mat[1][3] * det3_235_124 - mat[1][4] * det3_235_123;
	float det4_1235_1235 = mat[1][1] * det3_235_235 - mat[1][2] * det3_235_135 + mat[1][3] * det3_235_125 - mat[1][5] * det3_235_123;
	float det4_1235_1245 = mat[1][1] * det3_235_245 - mat[1][2] * det3_235_145 + mat[1][4] * det3_235_125 - mat[1][5] * det3_235_124;
	float det4_1235_1345 = mat[1][1] * det3_235_345 - mat[1][3] * det3_235_145 + mat[1][4] * det3_235_135 - mat[1][5] * det3_235_134;
	float det4_1235_2345 = mat[1][2] * det3_235_345 - mat[1][3] * det3_235_245 + mat[1][4] * det3_235_235 - mat[1][5] * det3_235_234;
	float det4_1245_0123 = mat[1][0] * det3_245_123 - mat[1][1] * det3_245_023 + mat[1][2] * det3_245_013 - mat[1][3] * det3_245_012;
	float det4_1245_0124 = mat[1][0] * det3_245_124 - mat[1][1] * det3_245_024 + mat[1][2] * det3_245_014 - mat[1][4] * det3_245_012;
	float det4_1245_0125 = mat[1][0] * det3_245_125 - mat[1][1] * det3_245_025 + mat[1][2] * det3_245_015 - mat[1][5] * det3_245_012;
	float det4_1245_0134 = mat[1][0] * det3_245_134 - mat[1][1] * det3_245_034 + mat[1][3] * det3_245_014 - mat[1][4] * det3_245_013;
	float det4_1245_0135 = mat[1][0] * det3_245_135 - mat[1][1] * det3_245_035 + mat[1][3] * det3_245_015 - mat[1][5] * det3_245_013;
	float det4_1245_0145 = mat[1][0] * det3_245_145 - mat[1][1] * det3_245_045 + mat[1][4] * det3_245_015 - mat[1][5] * det3_245_014;
	float det4_1245_0234 = mat[1][0] * det3_245_234 - mat[1][2] * det3_245_034 + mat[1][3] * det3_245_024 - mat[1][4] * det3_245_023;
	float det4_1245_0235 = mat[1][0] * det3_245_235 - mat[1][2] * det3_245_035 + mat[1][3] * det3_245_025 - mat[1][5] * det3_245_023;
	float det4_1245_0245 = mat[1][0] * det3_245_245 - mat[1][2] * det3_245_045 + mat[1][4] * det3_245_025 - mat[1][5] * det3_245_024;
	float det4_1245_0345 = mat[1][0] * det3_245_345 - mat[1][3] * det3_245_045 + mat[1][4] * det3_245_035 - mat[1][5] * det3_245_034;
	float det4_1245_1234 = mat[1][1] * det3_245_234 - mat[1][2] * det3_245_134 + mat[1][3] * det3_245_124 - mat[1][4] * det3_245_123;
	float det4_1245_1235 = mat[1][1] * det3_245_235 - mat[1][2] * det3_245_135 + mat[1][3] * det3_245_125 - mat[1][5] * det3_245_123;
	float det4_1245_1245 = mat[1][1] * det3_245_245 - mat[1][2] * det3_245_145 + mat[1][4] * det3_245_125 - mat[1][5] * det3_245_124;
	float det4_1245_1345 = mat[1][1] * det3_245_345 - mat[1][3] * det3_245_145 + mat[1][4] * det3_245_135 - mat[1][5] * det3_245_134;
	float det4_1245_2345 = mat[1][2] * det3_245_345 - mat[1][3] * det3_245_245 + mat[1][4] * det3_245_235 - mat[1][5] * det3_245_234;
	float det4_1345_0123 = mat[1][0] * det3_345_123 - mat[1][1] * det3_345_023 + mat[1][2] * det3_345_013 - mat[1][3] * det3_345_012;
	float det4_1345_0124 = mat[1][0] * det3_345_124 - mat[1][1] * det3_345_024 + mat[1][2] * det3_345_014 - mat[1][4] * det3_345_012;
	float det4_1345_0125 = mat[1][0] * det3_345_125 - mat[1][1] * det3_345_025 + mat[1][2] * det3_345_015 - mat[1][5] * det3_345_012;
	float det4_1345_0134 = mat[1][0] * det3_345_134 - mat[1][1] * det3_345_034 + mat[1][3] * det3_345_014 - mat[1][4] * det3_345_013;
	float det4_1345_0135 = mat[1][0] * det3_345_135 - mat[1][1] * det3_345_035 + mat[1][3] * det3_345_015 - mat[1][5] * det3_345_013;
	float det4_1345_0145 = mat[1][0] * det3_345_145 - mat[1][1] * det3_345_045 + mat[1][4] * det3_345_015 - mat[1][5] * det3_345_014;
	float det4_1345_0234 = mat[1][0] * det3_345_234 - mat[1][2] * det3_345_034 + mat[1][3] * det3_345_024 - mat[1][4] * det3_345_023;
	float det4_1345_0235 = mat[1][0] * det3_345_235 - mat[1][2] * det3_345_035 + mat[1][3] * det3_345_025 - mat[1][5] * det3_345_023;
	float det4_1345_0245 = mat[1][0] * det3_345_245 - mat[1][2] * det3_345_045 + mat[1][4] * det3_345_025 - mat[1][5] * det3_345_024;
	float det4_1345_0345 = mat[1][0] * det3_345_345 - mat[1][3] * det3_345_045 + mat[1][4] * det3_345_035 - mat[1][5] * det3_345_034;
	float det4_1345_1234 = mat[1][1] * det3_345_234 - mat[1][2] * det3_345_134 + mat[1][3] * det3_345_124 - mat[1][4] * det3_345_123;
	float det4_1345_1235 = mat[1][1] * det3_345_235 - mat[1][2] * det3_345_135 + mat[1][3] * det3_345_125 - mat[1][5] * det3_345_123;
	float det4_1345_1245 = mat[1][1] * det3_345_245 - mat[1][2] * det3_345_145 + mat[1][4] * det3_345_125 - mat[1][5] * det3_345_124;
	float det4_1345_1345 = mat[1][1] * det3_345_345 - mat[1][3] * det3_345_145 + mat[1][4] * det3_345_135 - mat[1][5] * det3_345_134;
	float det4_1345_2345 = mat[1][2] * det3_345_345 - mat[1][3] * det3_345_245 + mat[1][4] * det3_345_235 - mat[1][5] * det3_345_234;

	// remaining 5x5 sub-determinants
	float det5_01234_01234 = mat[0][0] * det4_1234_1234 - mat[0][1] * det4_1234_0234 + mat[0][2] * det4_1234_0134 - mat[0][3] * det4_1234_0124 + mat[0][4] * det4_1234_0123;
	float det5_01234_01235 = mat[0][0] * det4_1234_1235 - mat[0][1] * det4_1234_0235 + mat[0][2] * det4_1234_0135 - mat[0][3] * det4_1234_0125 + mat[0][5] * det4_1234_0123;
	float det5_01234_01245 = mat[0][0] * det4_1234_1245 - mat[0][1] * det4_1234_0245 + mat[0][2] * det4_1234_0145 - mat[0][4] * det4_1234_0125 + mat[0][5] * det4_1234_0124;
	float det5_01234_01345 = mat[0][0] * det4_1234_1345 - mat[0][1] * det4_1234_0345 + mat[0][3] * det4_1234_0145 - mat[0][4] * det4_1234_0135 + mat[0][5] * det4_1234_0134;
	float det5_01234_02345 = mat[0][0] * det4_1234_2345 - mat[0][2] * det4_1234_0345 + mat[0][3] * det4_1234_0245 - mat[0][4] * det4_1234_0235 + mat[0][5] * det4_1234_0234;
	float det5_01234_12345 = mat[0][1] * det4_1234_2345 - mat[0][2] * det4_1234_1345 + mat[0][3] * det4_1234_1245 - mat[0][4] * det4_1234_1235 + mat[0][5] * det4_1234_1234;
	float det5_01235_01234 = mat[0][0] * det4_1235_1234 - mat[0][1] * det4_1235_0234 + mat[0][2] * det4_1235_0134 - mat[0][3] * det4_1235_0124 + mat[0][4] * det4_1235_0123;
	float det5_01235_01235 = mat[0][0] * det4_1235_1235 - mat[0][1] * det4_1235_0235 + mat[0][2] * det4_1235_0135 - mat[0][3] * det4_1235_0125 + mat[0][5] * det4_1235_0123;
	float det5_01235_01245 = mat[0][0] * det4_1235_1245 - mat[0][1] * det4_1235_0245 + mat[0][2] * det4_1235_0145 - mat[0][4] * det4_1235_0125 + mat[0][5] * det4_1235_0124;
	float det5_01235_01345 = mat[0][0] * det4_1235_1345 - mat[0][1] * det4_1235_0345 + mat[0][3] * det4_1235_0145 - mat[0][4] * det4_1235_0135 + mat[0][5] * det4_1235_0134;
	float det5_01235_02345 = mat[0][0] * det4_1235_2345 - mat[0][2] * det4_1235_0345 + mat[0][3] * det4_1235_0245 - mat[0][4] * det4_1235_0235 + mat[0][5] * det4_1235_0234;
	float det5_01235_12345 = mat[0][1] * det4_1235_2345 - mat[0][2] * det4_1235_1345 + mat[0][3] * det4_1235_1245 - mat[0][4] * det4_1235_1235 + mat[0][5] * det4_1235_1234;
	float det5_01245_01234 = mat[0][0] * det4_1245_1234 - mat[0][1] * det4_1245_0234 + mat[0][2] * det4_1245_0134 - mat[0][3] * det4_1245_0124 + mat[0][4] * det4_1245_0123;
	float det5_01245_01235 = mat[0][0] * det4_1245_1235 - mat[0][1] * det4_1245_0235 + mat[0][2] * det4_1245_0135 - mat[0][3] * det4_1245_0125 + mat[0][5] * det4_1245_0123;
	float det5_01245_01245 = mat[0][0] * det4_1245_1245 - mat[0][1] * det4_1245_0245 + mat[0][2] * det4_1245_0145 - mat[0][4] * det4_1245_0125 + mat[0][5] * det4_1245_0124;
	float det5_01245_01345 = mat[0][0] * det4_1245_1345 - mat[0][1] * det4_1245_0345 + mat[0][3] * det4_1245_0145 - mat[0][4] * det4_1245_0135 + mat[0][5] * det4_1245_0134;
	float det5_01245_02345 = mat[0][0] * det4_1245_2345 - mat[0][2] * det4_1245_0345 + mat[0][3] * det4_1245_0245 - mat[0][4] * det4_1245_0235 + mat[0][5] * det4_1245_0234;
	float det5_01245_12345 = mat[0][1] * det4_1245_2345 - mat[0][2] * det4_1245_1345 + mat[0][3] * det4_1245_1245 - mat[0][4] * det4_1245_1235 + mat[0][5] * det4_1245_1234;
	float det5_01345_01234 = mat[0][0] * det4_1345_1234 - mat[0][1] * det4_1345_0234 + mat[0][2] * det4_1345_0134 - mat[0][3] * det4_1345_0124 + mat[0][4] * det4_1345_0123;
	float det5_01345_01235 = mat[0][0] * det4_1345_1235 - mat[0][1] * det4_1345_0235 + mat[0][2] * det4_1345_0135 - mat[0][3] * det4_1345_0125 + mat[0][5] * det4_1345_0123;
	float det5_01345_01245 = mat[0][0] * det4_1345_1245 - mat[0][1] * det4_1345_0245 + mat[0][2] * det4_1345_0145 - mat[0][4] * det4_1345_0125 + mat[0][5] * det4_1345_0124;
	float det5_01345_01345 = mat[0][0] * det4_1345_1345 - mat[0][1] * det4_1345_0345 + mat[0][3] * det4_1345_0145 - mat[0][4] * det4_1345_0135 + mat[0][5] * det4_1345_0134;
	float det5_01345_02345 = mat[0][0] * det4_1345_2345 - mat[0][2] * det4_1345_0345 + mat[0][3] * det4_1345_0245 - mat[0][4] * det4_1345_0235 + mat[0][5] * det4_1345_0234;
	float det5_01345_12345 = mat[0][1] * det4_1345_2345 - mat[0][2] * det4_1345_1345 + mat[0][3] * det4_1345_1245 - mat[0][4] * det4_1345_1235 + mat[0][5] * det4_1345_1234;
	float det5_02345_01234 = mat[0][0] * det4_2345_1234 - mat[0][1] * det4_2345_0234 + mat[0][2] * det4_2345_0134 - mat[0][3] * det4_2345_0124 + mat[0][4] * det4_2345_0123;
	float det5_02345_01235 = mat[0][0] * det4_2345_1235 - mat[0][1] * det4_2345_0235 + mat[0][2] * det4_2345_0135 - mat[0][3] * det4_2345_0125 + mat[0][5] * det4_2345_0123;
	float det5_02345_01245 = mat[0][0] * det4_2345_1245 - mat[0][1] * det4_2345_0245 + mat[0][2] * det4_2345_0145 - mat[0][4] * det4_2345_0125 + mat[0][5] * det4_2345_0124;
	float det5_02345_01345 = mat[0][0] * det4_2345_1345 - mat[0][1] * det4_2345_0345 + mat[0][3] * det4_2345_0145 - mat[0][4] * det4_2345_0135 + mat[0][5] * det4_2345_0134;
	float det5_02345_02345 = mat[0][0] * det4_2345_2345 - mat[0][2] * det4_2345_0345 + mat[0][3] * det4_2345_0245 - mat[0][4] * det4_2345_0235 + mat[0][5] * det4_2345_0234;
	float det5_02345_12345 = mat[0][1] * det4_2345_2345 - mat[0][2] * det4_2345_1345 + mat[0][3] * det4_2345_1245 - mat[0][4] * det4_2345_1235 + mat[0][5] * det4_2345_1234;

	mat[0][0] =  det5_12345_12345 * invDet;
	mat[0][1] = -det5_02345_12345 * invDet;
	mat[0][2] =  det5_01345_12345 * invDet;
	mat[0][3] = -det5_01245_12345 * invDet;
	mat[0][4] =  det5_01235_12345 * invDet;
	mat[0][5] = -det5_01234_12345 * invDet;

	mat[1][0] = -det5_12345_02345 * invDet;
	mat[1][1] =  det5_02345_02345 * invDet;
	mat[1][2] = -det5_01345_02345 * invDet;
	mat[1][3] =  det5_01245_02345 * invDet;
	mat[1][4] = -det5_01235_02345 * invDet;
	mat[1][5] =  det5_01234_02345 * invDet;

	mat[2][0] =  det5_12345_01345 * invDet;
	mat[2][1] = -det5_02345_01345 * invDet;
	mat[2][2] =  det5_01345_01345 * invDet;
	mat[2][3] = -det5_01245_01345 * invDet;
	mat[2][4] =  det5_01235_01345 * invDet;
	mat[2][5] = -det5_01234_01345 * invDet;

	mat[3][0] = -det5_12345_01245 * invDet;
	mat[3][1] =  det5_02345_01245 * invDet;
	mat[3][2] = -det5_01345_01245 * invDet;
	mat[3][3] =  det5_01245_01245 * invDet;
	mat[3][4] = -det5_01235_01245 * invDet;
	mat[3][5] =  det5_01234_01245 * invDet;

	mat[4][0] =  det5_12345_01235 * invDet;
	mat[4][1] = -det5_02345_01235 * invDet;
	mat[4][2] =  det5_01345_01235 * invDet;
	mat[4][3] = -det5_01245_01235 * invDet;
	mat[4][4] =  det5_01235_01235 * invDet;
	mat[4][5] = -det5_01234_01235 * invDet;

	mat[5][0] = -det5_12345_01234 * invDet;
	mat[5][1] =  det5_02345_01234 * invDet;
	mat[5][2] = -det5_01345_01234 * invDet;
	mat[5][3] =  det5_01245_01234 * invDet;
	mat[5][4] = -det5_01235_01234 * invDet;
	mat[5][5] =  det5_01234_01234 * invDet;

	return true;
#elif 0
	// 6*40 = 240 multiplications
	//			6 divisions
	float *mat = reinterpret_cast<float *>(this);
	float s;
	double d, di;

	di = mat[0];
	s = di;
	mat[0] = d = 1.0f / di;
	mat[1] *= d;
	mat[2] *= d;
	mat[3] *= d;
	mat[4] *= d;
	mat[5] *= d;
	d = -d;
	mat[6] *= d;
	mat[12] *= d;
	mat[18] *= d;
	mat[24] *= d;
	mat[30] *= d;
	d = mat[6] * di;
	mat[7] += mat[1] * d;
	mat[8] += mat[2] * d;
	mat[9] += mat[3] * d;
	mat[10] += mat[4] * d;
	mat[11] += mat[5] * d;
	d = mat[12] * di;
	mat[13] += mat[1] * d;
	mat[14] += mat[2] * d;
	mat[15] += mat[3] * d;
	mat[16] += mat[4] * d;
	mat[17] += mat[5] * d;
	d = mat[18] * di;
	mat[19] += mat[1] * d;
	mat[20] += mat[2] * d;
	mat[21] += mat[3] * d;
	mat[22] += mat[4] * d;
	mat[23] += mat[5] * d;
	d = mat[24] * di;
	mat[25] += mat[1] * d;
	mat[26] += mat[2] * d;
	mat[27] += mat[3] * d;
	mat[28] += mat[4] * d;
	mat[29] += mat[5] * d;
	d = mat[30] * di;
	mat[31] += mat[1] * d;
	mat[32] += mat[2] * d;
	mat[33] += mat[3] * d;
	mat[34] += mat[4] * d;
	mat[35] += mat[5] * d;
	di = mat[7];
	s *= di;
	mat[7] = d = 1.0f / di;
	mat[6] *= d;
	mat[8] *= d;
	mat[9] *= d;
	mat[10] *= d;
	mat[11] *= d;
	d = -d;
	mat[1] *= d;
	mat[13] *= d;
	mat[19] *= d;
	mat[25] *= d;
	mat[31] *= d;
	d = mat[1] * di;
	mat[0] += mat[6] * d;
	mat[2] += mat[8] * d;
	mat[3] += mat[9] * d;
	mat[4] += mat[10] * d;
	mat[5] += mat[11] * d;
	d = mat[13] * di;
	mat[12] += mat[6] * d;
	mat[14] += mat[8] * d;
	mat[15] += mat[9] * d;
	mat[16] += mat[10] * d;
	mat[17] += mat[11] * d;
	d = mat[19] * di;
	mat[18] += mat[6] * d;
	mat[20] += mat[8] * d;
	mat[21] += mat[9] * d;
	mat[22] += mat[10] * d;
	mat[23] += mat[11] * d;
	d = mat[25] * di;
	mat[24] += mat[6] * d;
	mat[26] += mat[8] * d;
	mat[27] += mat[9] * d;
	mat[28] += mat[10] * d;
	mat[29] += mat[11] * d;
	d = mat[31] * di;
	mat[30] += mat[6] * d;
	mat[32] += mat[8] * d;
	mat[33] += mat[9] * d;
	mat[34] += mat[10] * d;
	mat[35] += mat[11] * d;
	di = mat[14];
	s *= di;
	mat[14] = d = 1.0f / di;
	mat[12] *= d;
	mat[13] *= d;
	mat[15] *= d;
	mat[16] *= d;
	mat[17] *= d;
	d = -d;
	mat[2] *= d;
	mat[8] *= d;
	mat[20] *= d;
	mat[26] *= d;
	mat[32] *= d;
	d = mat[2] * di;
	mat[0] += mat[12] * d;
	mat[1] += mat[13] * d;
	mat[3] += mat[15] * d;
	mat[4] += mat[16] * d;
	mat[5] += mat[17] * d;
	d = mat[8] * di;
	mat[6] += mat[12] * d;
	mat[7] += mat[13] * d;
	mat[9] += mat[15] * d;
	mat[10] += mat[16] * d;
	mat[11] += mat[17] * d;
	d = mat[20] * di;
	mat[18] += mat[12] * d;
	mat[19] += mat[13] * d;
	mat[21] += mat[15] * d;
	mat[22] += mat[16] * d;
	mat[23] += mat[17] * d;
	d = mat[26] * di;
	mat[24] += mat[12] * d;
	mat[25] += mat[13] * d;
	mat[27] += mat[15] * d;
	mat[28] += mat[16] * d;
	mat[29] += mat[17] * d;
	d = mat[32] * di;
	mat[30] += mat[12] * d;
	mat[31] += mat[13] * d;
	mat[33] += mat[15] * d;
	mat[34] += mat[16] * d;
	mat[35] += mat[17] * d;
	di = mat[21];
	s *= di;
	mat[21] = d = 1.0f / di;
	mat[18] *= d;
	mat[19] *= d;
	mat[20] *= d;
	mat[22] *= d;
	mat[23] *= d;
	d = -d;
	mat[3] *= d;
	mat[9] *= d;
	mat[15] *= d;
	mat[27] *= d;
	mat[33] *= d;
	d = mat[3] * di;
	mat[0] += mat[18] * d;
	mat[1] += mat[19] * d;
	mat[2] += mat[20] * d;
	mat[4] += mat[22] * d;
	mat[5] += mat[23] * d;
	d = mat[9] * di;
	mat[6] += mat[18] * d;
	mat[7] += mat[19] * d;
	mat[8] += mat[20] * d;
	mat[10] += mat[22] * d;
	mat[11] += mat[23] * d;
	d = mat[15] * di;
	mat[12] += mat[18] * d;
	mat[13] += mat[19] * d;
	mat[14] += mat[20] * d;
	mat[16] += mat[22] * d;
	mat[17] += mat[23] * d;
	d = mat[27] * di;
	mat[24] += mat[18] * d;
	mat[25] += mat[19] * d;
	mat[26] += mat[20] * d;
	mat[28] += mat[22] * d;
	mat[29] += mat[23] * d;
	d = mat[33] * di;
	mat[30] += mat[18] * d;
	mat[31] += mat[19] * d;
	mat[32] += mat[20] * d;
	mat[34] += mat[22] * d;
	mat[35] += mat[23] * d;
	di = mat[28];
	s *= di;
	mat[28] = d = 1.0f / di;
	mat[24] *= d;
	mat[25] *= d;
	mat[26] *= d;
	mat[27] *= d;
	mat[29] *= d;
	d = -d;
	mat[4] *= d;
	mat[10] *= d;
	mat[16] *= d;
	mat[22] *= d;
	mat[34] *= d;
	d = mat[4] * di;
	mat[0] += mat[24] * d;
	mat[1] += mat[25] * d;
	mat[2] += mat[26] * d;
	mat[3] += mat[27] * d;
	mat[5] += mat[29] * d;
	d = mat[10] * di;
	mat[6] += mat[24] * d;
	mat[7] += mat[25] * d;
	mat[8] += mat[26] * d;
	mat[9] += mat[27] * d;
	mat[11] += mat[29] * d;
	d = mat[16] * di;
	mat[12] += mat[24] * d;
	mat[13] += mat[25] * d;
	mat[14] += mat[26] * d;
	mat[15] += mat[27] * d;
	mat[17] += mat[29] * d;
	d = mat[22] * di;
	mat[18] += mat[24] * d;
	mat[19] += mat[25] * d;
	mat[20] += mat[26] * d;
	mat[21] += mat[27] * d;
	mat[23] += mat[29] * d;
	d = mat[34] * di;
	mat[30] += mat[24] * d;
	mat[31] += mat[25] * d;
	mat[32] += mat[26] * d;
	mat[33] += mat[27] * d;
	mat[35] += mat[29] * d;
	di = mat[35];
	s *= di;
	mat[35] = d = 1.0f / di;
	mat[30] *= d;
	mat[31] *= d;
	mat[32] *= d;
	mat[33] *= d;
	mat[34] *= d;
	d = -d;
	mat[5] *= d;
	mat[11] *= d;
	mat[17] *= d;
	mat[23] *= d;
	mat[29] *= d;
	d = mat[5] * di;
	mat[0] += mat[30] * d;
	mat[1] += mat[31] * d;
	mat[2] += mat[32] * d;
	mat[3] += mat[33] * d;
	mat[4] += mat[34] * d;
	d = mat[11] * di;
	mat[6] += mat[30] * d;
	mat[7] += mat[31] * d;
	mat[8] += mat[32] * d;
	mat[9] += mat[33] * d;
	mat[10] += mat[34] * d;
	d = mat[17] * di;
	mat[12] += mat[30] * d;
	mat[13] += mat[31] * d;
	mat[14] += mat[32] * d;
	mat[15] += mat[33] * d;
	mat[16] += mat[34] * d;
	d = mat[23] * di;
	mat[18] += mat[30] * d;
	mat[19] += mat[31] * d;
	mat[20] += mat[32] * d;
	mat[21] += mat[33] * d;
	mat[22] += mat[34] * d;
	d = mat[29] * di;
	mat[24] += mat[30] * d;
	mat[25] += mat[31] * d;
	mat[26] += mat[32] * d;
	mat[27] += mat[33] * d;
	mat[28] += mat[34] * d;

	return ( s != 0.0f && !FLOAT_IS_NAN( s ) );
#else
	// 6*27+2*30 = 222 multiplications
	//		2*1  =	 2 divisions
	CMat3D r0, r1, r2, r3;
	float c0, c1, c2, det, invDet;
	float *mat = reinterpret_cast<float *>(this);

	// r0 = m0.Inverse();
	c0 = mat[1*6+1] * mat[2*6+2] - mat[1*6+2] * mat[2*6+1];
	c1 = mat[1*6+2] * mat[2*6+0] - mat[1*6+0] * mat[2*6+2];
	c2 = mat[1*6+0] * mat[2*6+1] - mat[1*6+1] * mat[2*6+0];

	det = mat[0*6+0] * c0 + mat[0*6+1] * c1 + mat[0*6+2] * c2;

	if ( CMath::fabs( det ) < MATRIX_INVERSE_EPSILON ) {
		return false;
	}

	invDet = 1.0f / det;

	r0[0][0] = c0 * invDet;
	r0[0][1] = ( mat[0*6+2] * mat[2*6+1] - mat[0*6+1] * mat[2*6+2] ) * invDet;
	r0[0][2] = ( mat[0*6+1] * mat[1*6+2] - mat[0*6+2] * mat[1*6+1] ) * invDet;
	r0[1][0] = c1 * invDet;
	r0[1][1] = ( mat[0*6+0] * mat[2*6+2] - mat[0*6+2] * mat[2*6+0] ) * invDet;
	r0[1][2] = ( mat[0*6+2] * mat[1*6+0] - mat[0*6+0] * mat[1*6+2] ) * invDet;
	r0[2][0] = c2 * invDet;
	r0[2][1] = ( mat[0*6+1] * mat[2*6+0] - mat[0*6+0] * mat[2*6+1] ) * invDet;
	r0[2][2] = ( mat[0*6+0] * mat[1*6+1] - mat[0*6+1] * mat[1*6+0] ) * invDet;

	// r1 = r0 * m1;
	r1[0][0] = r0[0][0] * mat[0*6+3] + r0[0][1] * mat[1*6+3] + r0[0][2] * mat[2*6+3];
	r1[0][1] = r0[0][0] * mat[0*6+4] + r0[0][1] * mat[1*6+4] + r0[0][2] * mat[2*6+4];
	r1[0][2] = r0[0][0] * mat[0*6+5] + r0[0][1] * mat[1*6+5] + r0[0][2] * mat[2*6+5];
	r1[1][0] = r0[1][0] * mat[0*6+3] + r0[1][1] * mat[1*6+3] + r0[1][2] * mat[2*6+3];
	r1[1][1] = r0[1][0] * mat[0*6+4] + r0[1][1] * mat[1*6+4] + r0[1][2] * mat[2*6+4];
	r1[1][2] = r0[1][0] * mat[0*6+5] + r0[1][1] * mat[1*6+5] + r0[1][2] * mat[2*6+5];
	r1[2][0] = r0[2][0] * mat[0*6+3] + r0[2][1] * mat[1*6+3] + r0[2][2] * mat[2*6+3];
	r1[2][1] = r0[2][0] * mat[0*6+4] + r0[2][1] * mat[1*6+4] + r0[2][2] * mat[2*6+4];
	r1[2][2] = r0[2][0] * mat[0*6+5] + r0[2][1] * mat[1*6+5] + r0[2][2] * mat[2*6+5];

	// r2 = m2 * r1;
	r2[0][0] = mat[3*6+0] * r1[0][0] + mat[3*6+1] * r1[1][0] + mat[3*6+2] * r1[2][0];
	r2[0][1] = mat[3*6+0] * r1[0][1] + mat[3*6+1] * r1[1][1] + mat[3*6+2] * r1[2][1];
	r2[0][2] = mat[3*6+0] * r1[0][2] + mat[3*6+1] * r1[1][2] + mat[3*6+2] * r1[2][2];
	r2[1][0] = mat[4*6+0] * r1[0][0] + mat[4*6+1] * r1[1][0] + mat[4*6+2] * r1[2][0];
	r2[1][1] = mat[4*6+0] * r1[0][1] + mat[4*6+1] * r1[1][1] + mat[4*6+2] * r1[2][1];
	r2[1][2] = mat[4*6+0] * r1[0][2] + mat[4*6+1] * r1[1][2] + mat[4*6+2] * r1[2][2];
	r2[2][0] = mat[5*6+0] * r1[0][0] + mat[5*6+1] * r1[1][0] + mat[5*6+2] * r1[2][0];
	r2[2][1] = mat[5*6+0] * r1[0][1] + mat[5*6+1] * r1[1][1] + mat[5*6+2] * r1[2][1];
	r2[2][2] = mat[5*6+0] * r1[0][2] + mat[5*6+1] * r1[1][2] + mat[5*6+2] * r1[2][2];

	// r3 = r2 - m3;
	r3[0][0] = r2[0][0] - mat[3*6+3];
	r3[0][1] = r2[0][1] - mat[3*6+4];
	r3[0][2] = r2[0][2] - mat[3*6+5];
	r3[1][0] = r2[1][0] - mat[4*6+3];
	r3[1][1] = r2[1][1] - mat[4*6+4];
	r3[1][2] = r2[1][2] - mat[4*6+5];
	r3[2][0] = r2[2][0] - mat[5*6+3];
	r3[2][1] = r2[2][1] - mat[5*6+4];
	r3[2][2] = r2[2][2] - mat[5*6+5];

	// r3.InverseSelf();
	r2[0][0] = r3[1][1] * r3[2][2] - r3[1][2] * r3[2][1];
	r2[1][0] = r3[1][2] * r3[2][0] - r3[1][0] * r3[2][2];
	r2[2][0] = r3[1][0] * r3[2][1] - r3[1][1] * r3[2][0];

	det = r3[0][0] * r2[0][0] + r3[0][1] * r2[1][0] + r3[0][2] * r2[2][0];

	if ( CMath::fabs( det ) < MATRIX_INVERSE_EPSILON ) {
		return false;
	}

	invDet = 1.0f / det;

	r2[0][1] = r3[0][2] * r3[2][1] - r3[0][1] * r3[2][2];
	r2[0][2] = r3[0][1] * r3[1][2] - r3[0][2] * r3[1][1];
	r2[1][1] = r3[0][0] * r3[2][2] - r3[0][2] * r3[2][0];
	r2[1][2] = r3[0][2] * r3[1][0] - r3[0][0] * r3[1][2];
	r2[2][1] = r3[0][1] * r3[2][0] - r3[0][0] * r3[2][1];
	r2[2][2] = r3[0][0] * r3[1][1] - r3[0][1] * r3[1][0];

	r3[0][0] = r2[0][0] * invDet;
	r3[0][1] = r2[0][1] * invDet;
	r3[0][2] = r2[0][2] * invDet;
	r3[1][0] = r2[1][0] * invDet;
	r3[1][1] = r2[1][1] * invDet;
	r3[1][2] = r2[1][2] * invDet;
	r3[2][0] = r2[2][0] * invDet;
	r3[2][1] = r2[2][1] * invDet;
	r3[2][2] = r2[2][2] * invDet;

	// r2 = m2 * r0;
	r2[0][0] = mat[3*6+0] * r0[0][0] + mat[3*6+1] * r0[1][0] + mat[3*6+2] * r0[2][0];
	r2[0][1] = mat[3*6+0] * r0[0][1] + mat[3*6+1] * r0[1][1] + mat[3*6+2] * r0[2][1];
	r2[0][2] = mat[3*6+0] * r0[0][2] + mat[3*6+1] * r0[1][2] + mat[3*6+2] * r0[2][2];
	r2[1][0] = mat[4*6+0] * r0[0][0] + mat[4*6+1] * r0[1][0] + mat[4*6+2] * r0[2][0];
	r2[1][1] = mat[4*6+0] * r0[0][1] + mat[4*6+1] * r0[1][1] + mat[4*6+2] * r0[2][1];
	r2[1][2] = mat[4*6+0] * r0[0][2] + mat[4*6+1] * r0[1][2] + mat[4*6+2] * r0[2][2];
	r2[2][0] = mat[5*6+0] * r0[0][0] + mat[5*6+1] * r0[1][0] + mat[5*6+2] * r0[2][0];
	r2[2][1] = mat[5*6+0] * r0[0][1] + mat[5*6+1] * r0[1][1] + mat[5*6+2] * r0[2][1];
	r2[2][2] = mat[5*6+0] * r0[0][2] + mat[5*6+1] * r0[1][2] + mat[5*6+2] * r0[2][2];

	// m2 = r3 * r2;
	mat[3*6+0] = r3[0][0] * r2[0][0] + r3[0][1] * r2[1][0] + r3[0][2] * r2[2][0];
	mat[3*6+1] = r3[0][0] * r2[0][1] + r3[0][1] * r2[1][1] + r3[0][2] * r2[2][1];
	mat[3*6+2] = r3[0][0] * r2[0][2] + r3[0][1] * r2[1][2] + r3[0][2] * r2[2][2];
	mat[4*6+0] = r3[1][0] * r2[0][0] + r3[1][1] * r2[1][0] + r3[1][2] * r2[2][0];
	mat[4*6+1] = r3[1][0] * r2[0][1] + r3[1][1] * r2[1][1] + r3[1][2] * r2[2][1];
	mat[4*6+2] = r3[1][0] * r2[0][2] + r3[1][1] * r2[1][2] + r3[1][2] * r2[2][2];
	mat[5*6+0] = r3[2][0] * r2[0][0] + r3[2][1] * r2[1][0] + r3[2][2] * r2[2][0];
	mat[5*6+1] = r3[2][0] * r2[0][1] + r3[2][1] * r2[1][1] + r3[2][2] * r2[2][1];
	mat[5*6+2] = r3[2][0] * r2[0][2] + r3[2][1] * r2[1][2] + r3[2][2] * r2[2][2];

	// m0 = r0 - r1 * m2;
	mat[0*6+0] = r0[0][0] - r1[0][0] * mat[3*6+0] - r1[0][1] * mat[4*6+0] - r1[0][2] * mat[5*6+0];
	mat[0*6+1] = r0[0][1] - r1[0][0] * mat[3*6+1] - r1[0][1] * mat[4*6+1] - r1[0][2] * mat[5*6+1];
	mat[0*6+2] = r0[0][2] - r1[0][0] * mat[3*6+2] - r1[0][1] * mat[4*6+2] - r1[0][2] * mat[5*6+2];
	mat[1*6+0] = r0[1][0] - r1[1][0] * mat[3*6+0] - r1[1][1] * mat[4*6+0] - r1[1][2] * mat[5*6+0];
	mat[1*6+1] = r0[1][1] - r1[1][0] * mat[3*6+1] - r1[1][1] * mat[4*6+1] - r1[1][2] * mat[5*6+1];
	mat[1*6+2] = r0[1][2] - r1[1][0] * mat[3*6+2] - r1[1][1] * mat[4*6+2] - r1[1][2] * mat[5*6+2];
	mat[2*6+0] = r0[2][0] - r1[2][0] * mat[3*6+0] - r1[2][1] * mat[4*6+0] - r1[2][2] * mat[5*6+0];
	mat[2*6+1] = r0[2][1] - r1[2][0] * mat[3*6+1] - r1[2][1] * mat[4*6+1] - r1[2][2] * mat[5*6+1];
	mat[2*6+2] = r0[2][2] - r1[2][0] * mat[3*6+2] - r1[2][1] * mat[4*6+2] - r1[2][2] * mat[5*6+2];

	// m1 = r1 * r3;
	mat[0*6+3] = r1[0][0] * r3[0][0] + r1[0][1] * r3[1][0] + r1[0][2] * r3[2][0];
	mat[0*6+4] = r1[0][0] * r3[0][1] + r1[0][1] * r3[1][1] + r1[0][2] * r3[2][1];
	mat[0*6+5] = r1[0][0] * r3[0][2] + r1[0][1] * r3[1][2] + r1[0][2] * r3[2][2];
	mat[1*6+3] = r1[1][0] * r3[0][0] + r1[1][1] * r3[1][0] + r1[1][2] * r3[2][0];
	mat[1*6+4] = r1[1][0] * r3[0][1] + r1[1][1] * r3[1][1] + r1[1][2] * r3[2][1];
	mat[1*6+5] = r1[1][0] * r3[0][2] + r1[1][1] * r3[1][2] + r1[1][2] * r3[2][2];
	mat[2*6+3] = r1[2][0] * r3[0][0] + r1[2][1] * r3[1][0] + r1[2][2] * r3[2][0];
	mat[2*6+4] = r1[2][0] * r3[0][1] + r1[2][1] * r3[1][1] + r1[2][2] * r3[2][1];
	mat[2*6+5] = r1[2][0] * r3[0][2] + r1[2][1] * r3[1][2] + r1[2][2] * r3[2][2];

	// m3 = -r3;
	mat[3*6+3] = -r3[0][0];
	mat[3*6+4] = -r3[0][1];
	mat[3*6+5] = -r3[0][2];
	mat[4*6+3] = -r3[1][0];
	mat[4*6+4] = -r3[1][1];
	mat[4*6+5] = -r3[1][2];
	mat[5*6+3] = -r3[2][0];
	mat[5*6+4] = -r3[2][1];
	mat[5*6+5] = -r3[2][2];

	return true;
#endif
}

/*
=============
CMat6D::toString
=============
*/
const char *CMat6D::toString( int precision ) const {
	return CMyString::floatArraytoString( toFloatPtr(), getDimension(), precision );
}


//===============================================================
//
//  CMatXD
//
//===============================================================

float	CMatXD::temp[MATX_MAX_TEMP+4];
float *	CMatXD::tempPtr = (float *) ( ( (uintptr_t) CMatXD::temp + 15 ) & ~15 );
int		CMatXD::tempIndex = 0;


/*
============
CMatXD::changeSize
============
*/
void CMatXD::changeSize( int rows, int columns, bool makeZero ) {
	int alloc = ( rows * columns + 3 ) & ~3;
	if ( alloc > alloced && alloced != -1 ) {
		float *oldMat = mat;
		mat = (float *) mem_Alloc16( alloc * sizeof( float ) );
		if ( makeZero ) {
			memset( mat, 0, alloc * sizeof( float ) );
		}
		alloced = alloc;
		if ( oldMat ) {
			int minRow = MIN( numRows, rows );
			int minColumn = MIN( numColumns, columns );
			for ( int i = 0; i < minRow; i++ ) {
				for ( int j = 0; j < minColumn; j++ ) {
					mat[ i * columns + j ] = oldMat[ i * numColumns + j ];
				}
			}
			mem_Free16( oldMat );
		}
	} else {
		if ( columns < numColumns ) {
			int minRow = MIN( numRows, rows );
			for ( int i = 0; i < minRow; i++ ) {
				for ( int j = 0; j < columns; j++ ) {
					mat[ i * columns + j ] = mat[ i * numColumns + j ];
				}
			}
		} else if ( columns > numColumns ) {
			for ( int i = MIN( numRows, rows ) - 1; i >= 0; i-- ) {
				if ( makeZero ) {
					for ( int j = columns - 1; j >= numColumns; j-- ) {
						mat[ i * columns + j ] = 0.0f;
					}
				}
				for ( int j = numColumns - 1; j >= 0; j-- ) {
					mat[ i * columns + j ] = mat[ i * numColumns + j ];
				}
			}
		}
		if ( makeZero && rows > numRows ) {
			memset( mat + numRows * columns, 0, ( rows - numRows ) * columns * sizeof( float ) );
		}
	}
	numRows = rows;
	numColumns = columns;
	MATX_CLEAREND();
}

/*
============
CMatXD::removeRow
============
*/
CMatXD &CMatXD::removeRow( int r ) {
	int i;

	SMF_ASSERT( r < numRows );

	numRows--;

	for ( i = r; i < numRows; i++ ) {
		memcpy( &mat[i * numColumns], &mat[( i + 1 ) * numColumns], numColumns * sizeof( float ) );
	}

	return *this;
}

/*
============
CMatXD::removeColumn
============
*/
CMatXD &CMatXD::removeColumn( int r ) {
	int i;

	SMF_ASSERT( r < numColumns );

	numColumns--;

	for ( i = 0; i < numRows - 1; i++ ) {
		memmove( &mat[i * numColumns + r], &mat[i * ( numColumns + 1 ) + r + 1], numColumns * sizeof( float ) );
	}
	memmove( &mat[i * numColumns + r], &mat[i * ( numColumns + 1 ) + r + 1], ( numColumns - r ) * sizeof( float ) );

	return *this;
}

/*
============
CMatXD::removeRowColumn
============
*/
CMatXD &CMatXD::removeRowColumn( int r ) {
	int i;

	SMF_ASSERT( r < numRows && r < numColumns );

	numRows--;
	numColumns--;

	if ( r > 0 ) {
		for ( i = 0; i < r - 1; i++ ) {
			memmove( &mat[i * numColumns + r], &mat[i * ( numColumns + 1 ) + r + 1], numColumns * sizeof( float ) );
		}
		memmove( &mat[i * numColumns + r], &mat[i * ( numColumns + 1 ) + r + 1], ( numColumns - r ) * sizeof( float ) );
	}

	memcpy( &mat[r * numColumns], &mat[( r + 1 ) * ( numColumns + 1 )], r * sizeof( float ) );

	for ( i = r; i < numRows - 1; i++ ) {
		memcpy( &mat[i * numColumns + r], &mat[( i + 1 ) * ( numColumns + 1 ) + r + 1], numColumns * sizeof( float ) );
	}
	memcpy( &mat[i * numColumns + r], &mat[( i + 1 ) * ( numColumns + 1 ) + r + 1], ( numColumns - r ) * sizeof( float ) );

	return *this;
}

/*
============
CMatXD::isOrthogonal

  returns true if (*this) * this->transpose() == identity
============
*/
bool CMatXD::isOrthogonal( const float epsilon ) const {
	float *ptr1, *ptr2, sum;

	if ( !isSquare() ) {
		return false;
	}

	ptr1 = mat;
	for ( int i = 0; i < numRows; i++ ) {
		for ( int j = 0; j < numColumns; j++ ) {
			ptr2 = mat + j;
			sum = ptr1[0] * ptr2[0] - (float) ( i == j );
			for ( int n = 1; n < numColumns; n++ ) {
				ptr2 += numColumns;
				sum += ptr1[n] * ptr2[0];
			}
			if ( CMath::fabs( sum ) > epsilon ) {
				return false;
			}
		}
		ptr1 += numColumns;
	}
	return true;
}

/*
============
CMatXD::isOrthonormal

  returns true if (*this) * this->transpose() == identity and the length of each column vector is 1
============
*/
bool CMatXD::isOrthonormal( const float epsilon ) const {
	float *ptr1, *ptr2, sum;

	if ( !isSquare() ) {
		return false;
	}

	ptr1 = mat;
	for ( int i = 0; i < numRows; i++ ) {
		for ( int j = 0; j < numColumns; j++ ) {
			ptr2 = mat + j;
			sum = ptr1[0] * ptr2[0] - (float) ( i == j );
			for ( int n = 1; n < numColumns; n++ ) {
				ptr2 += numColumns;
				sum += ptr1[n] * ptr2[0];
			}
			if ( CMath::fabs( sum ) > epsilon ) {
				return false;
			}
		}
		ptr1 += numColumns;

		ptr2 = mat + i;
		sum = ptr2[0] * ptr2[0] - 1.0f;
		for ( i = 1; i < numRows; i++ ) {
			ptr2 += numColumns;
			sum += ptr2[i] * ptr2[i];
		}
		if ( CMath::fabs( sum ) > epsilon ) {
			return false;
		}
	}
	return true;
}

/*
============
CMatXD::isPMatrix

  returns true if the matrix is a P-matrix
  A square matrix is a P-matrix if all its principal minors are positive.
============
*/
bool CMatXD::isPMatrix( const float epsilon ) const {
	int i, j;
	float d;
	CMatXD m;

	if ( !isSquare() ) {
		return false;
	}

	if ( numRows <= 0 ) {
		return true;
	}

	if ( (*this)[0][0] <= epsilon ) {
		return false;
	}

	if ( numRows <= 1 ) {
		return true;
	}

	m.setData( numRows - 1, numColumns - 1, MATX_ALLOCAFLOAT( ( numRows - 1 ) * ( numColumns - 1 ) ) );

	for ( i = 1; i < numRows; i++ ) {
		for ( j = 1; j < numColumns; j++ ) {
			m[i-1][j-1] = (*this)[i][j];
		}
	}

	if ( !m.isPMatrix( epsilon ) ) {
		return false;
	}

	for ( i = 1; i < numRows; i++ ) {
		d = (*this)[i][0] / (*this)[0][0];
		for ( j = 1; j < numColumns; j++ ) {
			m[i-1][j-1] = (*this)[i][j] - d * (*this)[0][j];
		}
	}

	if ( !m.isPMatrix( epsilon ) ) {
		return false;
	}

	return true;
}

/*
============
CMatXD::isZMatrix

  returns true if the matrix is a Z-matrix
  A square matrix M is a Z-matrix if M[i][j] <= 0 for all i != j.
============
*/
bool CMatXD::isZMatrix( const float epsilon ) const {
	int i, j;

	if ( !isSquare() ) {
		return false;
	}

	for ( i = 0; i < numRows; i++ ) {
		for ( j = 0; j < numColumns; j++ ) {
			if ( (*this)[i][j] > epsilon && i != j ) {
				return false;
			}
		}
	}
	return true;
}

/*
============
CMatXD::isPositiveDefinite

  returns true if the matrix is Positive Definite (PD)
  A square matrix M of order n is said to be PD if y'My > 0 for all vectors y of dimension n, y != 0.
============
*/
bool CMatXD::isPositiveDefinite( const float epsilon ) const {
	int i, j, k;
	float d, s;
	CMatXD m;

	// the matrix must be square
	if ( !isSquare() ) {
		return false;
	}

	// copy matrix
	m.setData( numRows, numColumns, MATX_ALLOCAFLOAT( numRows * numColumns ) );
	m = *this;

	// add transpose
	for ( i = 0; i < numRows; i++ ) {
		for ( j = 0; j < numColumns; j++ ) {
			m[i][j] += (*this)[j][i];
		}
	}

	// test Positive Definiteness with Gaussian pivot steps
	for ( i = 0; i < numRows; i++ ) {

		for ( j = i; j < numColumns; j++ ) {
			if ( m[j][j] <= epsilon ) {
				return false;
			}
		}

		d = 1.0f / m[i][i];
		for ( j = i + 1; j < numColumns; j++ ) {
			s = d * m[j][i];
			m[j][i] = 0.0f;
			for ( k = i + 1; k < numRows; k++ ) {
				m[j][k] -= s * m[i][k];
			}
		}
	}

	return true;
}

/*
============
CMatXD::isSymmetricPositiveDefinite

  returns true if the matrix is Symmetric Positive Definite (PD)
============
*/
bool CMatXD::isSymmetricPositiveDefinite( const float epsilon ) const {
	CMatXD m;

	// the matrix must be symmetric
	if ( !isSymmetric( epsilon ) ) {
		return false;
	}

	// copy matrix
	m.setData( numRows, numColumns, MATX_ALLOCAFLOAT( numRows * numColumns ) );
	m = *this;

	// being able to obtain Cholesky factors is both a necessary and sufficient condition for positive definiteness
	return m.cholesky_Factor();
}

/*
============
CMatXD::isPositiveSemiDefinite

  returns true if the matrix is Positive Semi Definite (PSD)
  A square matrix M of order n is said to be PSD if y'My >= 0 for all vectors y of dimension n, y != 0.
============
*/
bool CMatXD::isPositiveSemiDefinite( const float epsilon ) const {
	int i, j, k;
	float d, s;
	CMatXD m;

	// the matrix must be square
	if ( !isSquare() ) {
		return false;
	}

	// copy original matrix
	m.setData( numRows, numColumns, MATX_ALLOCAFLOAT( numRows * numColumns ) );
	m = *this;

	// add transpose
	for ( i = 0; i < numRows; i++ ) {
		for ( j = 0; j < numColumns; j++ ) {
			m[i][j] += (*this)[j][i];
		}
	}

	// test Positive Semi Definiteness with Gaussian pivot steps
	for ( i = 0; i < numRows; i++ ) {

		for ( j = i; j < numColumns; j++ ) {
			if ( m[j][j] < -epsilon ) {
				return false;
			}
			if ( m[j][j] > epsilon ) {
				continue;
			}
			for ( k = 0; k < numRows; k++ ) {
				if ( CMath::fabs( m[k][j] ) > epsilon ) {
					return false;
				}
				if ( CMath::fabs( m[j][k] ) > epsilon ) {
					return false;
				}
			}
		}

		if ( m[i][i] <= epsilon ) {
			continue;
		}

		d = 1.0f / m[i][i];
		for ( j = i + 1; j < numColumns; j++ ) {
			s = d * m[j][i];
			m[j][i] = 0.0f;
			for ( k = i + 1; k < numRows; k++ ) {
				m[j][k] -= s * m[i][k];
			}
		}
	}

	return true;
}

/*
============
CMatXD::isSymmetricPositiveSemiDefinite

  returns true if the matrix is Symmetric Positive Semi Definite (PSD)
============
*/
bool CMatXD::isSymmetricPositiveSemiDefinite( const float epsilon ) const {

	// the matrix must be symmetric
	if ( !isSymmetric( epsilon ) ) {
		return false;
	}

	return isPositiveSemiDefinite( epsilon );
}

/*
============
CMatXD::lowerTriangularInverse

  in-place inversion of the lower triangular matrix
============
*/
bool CMatXD::lowerTriangularInverse() {
	int i, j, k;
	double d, sum;

	for ( i = 0; i < numRows; i++ ) {
		d = (*this)[i][i];
		if ( d == 0.0f ) {
			return false;
		}
		(*this)[i][i] = d = 1.0f / d;

		for ( j = 0; j < i; j++ ) {
			sum = 0.0f;
			for ( k = j; k < i; k++ ) {
				sum -= (*this)[i][k] * (*this)[k][j];
			}
			(*this)[i][j] = sum * d;
		}
	}
	return true;
}

/*
============
CMatXD::upperTriangularInverse

  in-place inversion of the upper triangular matrix
============
*/
bool CMatXD::upperTriangularInverse() {
	int i, j, k;
	double d, sum;

	for ( i = numRows-1; i >= 0; i-- ) {
		d = (*this)[i][i];
		if ( d == 0.0f ) {
			return false;
		}
		(*this)[i][i] = d = 1.0f / d;

		for ( j = numRows-1; j > i; j-- ) {
			sum = 0.0f;
			for ( k = j; k > i; k-- ) {
				sum -= (*this)[i][k] * (*this)[k][j];
			}
			(*this)[i][j] = sum * d;
		}
	}
	return true;
}

/*
=============
CMatXD::toString
=============
*/
const char *CMatXD::toString( int precision ) const {
	return CMyString::floatArraytoString( toFloatPtr(), getDimension(), precision );
}

/*
============
CMatXD::update_RankOne

  Updates the matrix to obtain the matrix: A + alpha * v * w'
============
*/
void CMatXD::update_RankOne( const CVecXD &v, const CVecXD &w, float alpha ) {
	int i, j;
	float s;

	SMF_ASSERT( v.getSize() >= numRows );
	SMF_ASSERT( w.getSize() >= numColumns );

	for ( i = 0; i < numRows; i++ ) {
		s = alpha * v[i];
		for ( j = 0; j < numColumns; j++ ) {
			(*this)[i][j] += s * w[j];
		}
	}
}

/*
============
CMatXD::update_RankOneSymmetric

  Updates the matrix to obtain the matrix: A + alpha * v * v'
============
*/
void CMatXD::update_RankOneSymmetric( const CVecXD &v, float alpha ) {
	int i, j;
	float s;

	SMF_ASSERT( numRows == numColumns );
	SMF_ASSERT( v.getSize() >= numRows );

	for ( i = 0; i < numRows; i++ ) {
		s = alpha * v[i];
		for ( j = 0; j < numColumns; j++ ) {
			(*this)[i][j] += s * v[j];
		}
	}
}

/*
============
CMatXD::update_RowColumn

  Updates the matrix to obtain the matrix:

      [ 0  a  0 ]
  A + [ d  b  e ]
      [ 0  c  0 ]

  where: a = v[0,r-1], b = v[r], c = v[r+1,numRows-1], d = w[0,r-1], w[r] = 0.0f, e = w[r+1,numColumns-1]
============
*/
void CMatXD::update_RowColumn( const CVecXD &v, const CVecXD &w, int r ) {
	int i;

	SMF_ASSERT( w[r] == 0.0f );
	SMF_ASSERT( v.getSize() >= numColumns );
	SMF_ASSERT( w.getSize() >= numRows );

	for ( i = 0; i < numRows; i++ ) {
		(*this)[i][r] += v[i];
	}
	for ( i = 0; i < numColumns; i++ ) {
		(*this)[r][i] += w[i];
	}
}

/*
============
CMatXD::update_RowColumnSymmetric

  Updates the matrix to obtain the matrix:

      [ 0  a  0 ]
  A + [ a  b  c ]
      [ 0  c  0 ]

  where: a = v[0,r-1], b = v[r], c = v[r+1,numRows-1]
============
*/
void CMatXD::update_RowColumnSymmetric( const CVecXD &v, int r ) {
	int i;

	SMF_ASSERT( numRows == numColumns );
	SMF_ASSERT( v.getSize() >= numRows );

	for ( i = 0; i < r; i++ ) {
		(*this)[i][r] += v[i];
		(*this)[r][i] += v[i];
	}
	(*this)[r][r] += v[r];
	for ( i = r+1; i < numRows; i++ ) {
		(*this)[i][r] += v[i];
		(*this)[r][i] += v[i];
	}
}

/*
============
CMatXD::update_Increment

  Updates the matrix to obtain the matrix:

  [ A  a ]
  [ c  b ]

  where: a = v[0,numRows-1], b = v[numRows], c = w[0,numColumns-1]], w[numColumns] = 0
============
*/
void CMatXD::update_Increment( const CVecXD &v, const CVecXD &w ) {
	int i;

	SMF_ASSERT( numRows == numColumns );
	SMF_ASSERT( v.getSize() >= numRows+1 );
	SMF_ASSERT( w.getSize() >= numColumns+1 );

	changeSize( numRows+1, numColumns+1, false );

	for ( i = 0; i < numRows; i++ ) {
		(*this)[i][numColumns-1] = v[i];
	}
	for ( i = 0; i < numColumns-1; i++ ) {
		(*this)[numRows-1][i] = w[i];
	}
}

/*
============
CMatXD::update_IncrementSymmetric

  Updates the matrix to obtain the matrix:

  [ A  a ]
  [ a  b ]

  where: a = v[0,numRows-1], b = v[numRows]
============
*/
void CMatXD::update_IncrementSymmetric( const CVecXD &v ) {
	int i;

	SMF_ASSERT( numRows == numColumns );
	SMF_ASSERT( v.getSize() >= numRows+1 );

	changeSize( numRows+1, numColumns+1, false );

	for ( i = 0; i < numRows-1; i++ ) {
		(*this)[i][numColumns-1] = v[i];
	}
	for ( i = 0; i < numColumns; i++ ) {
		(*this)[numRows-1][i] = v[i];
	}
}

/*
============
CMatXD::update_Decrement

  Updates the matrix to obtain a matrix with row r and column r removed.
============
*/
void CMatXD::update_Decrement( int r ) {
	removeRowColumn( r );
}

/*
============
CMatXD::inverse_GaussJordan

  in-place inversion using Gauss-Jordan elimination
============
*/
bool CMatXD::inverse_GaussJordan() {
	int i, j, k, r, c;
	float d, max;

	SMF_ASSERT( numRows == numColumns );

	int *columnIndex = (int *) _allocafloat16( numRows * sizeof( int ) );
	int *rowIndex = (int *) _allocafloat16( numRows * sizeof( int ) );
	bool *pivot = (bool *) _allocafloat16( numRows * sizeof( bool ) );

	memset( pivot, 0, numRows * sizeof( bool ) );

	// elimination with full pivoting
	for ( i = 0; i < numRows; i++ ) {

		// search the whole matrix except for pivoted rows for the maximum absolute value
		max = 0.0f;
		r = c = 0;
		for ( j = 0; j < numRows; j++ ) {
			if ( !pivot[j] ) {
				for ( k = 0; k < numRows; k++ ) {
					if ( !pivot[k] ) {
						d = CMath::fabs( (*this)[j][k] );
						if ( d > max ) {
							max = d;
							r = j;
							c = k;
						}
					}
				}
			}
		}

		if ( max == 0.0f ) {
			// matrix is not invertible
			return false;
		}

		pivot[c] = true;

		// swap rows such that entry (c,c) has the pivot entry
		if ( r != c ) {
			swapRows( r, c );
		}

		// keep track of the row permutation
		rowIndex[i] = r;
		columnIndex[i] = c;

		// scale the row to make the pivot entry equal to 1
		d = 1.0f / (*this)[c][c];
		(*this)[c][c] = 1.0f;
		for ( k = 0; k < numRows; k++ ) {
			(*this)[c][k] *= d;
		}

		// zero out the pivot column entries in the other rows
		for ( j = 0; j < numRows; j++ ) {
			if ( j != c ) {
				d = (*this)[j][c];
				(*this)[j][c] = 0.0f;
				for ( k = 0; k < numRows; k++ ) {
					(*this)[j][k] -= (*this)[c][k] * d;
				}
			}
		}
	}

	// reorder rows to store the inverse of the original matrix
	for ( j = numRows - 1; j >= 0; j-- ) {
		if ( rowIndex[j] != columnIndex[j] ) {
			for ( k = 0; k < numRows; k++ ) {
				d = (*this)[k][rowIndex[j]];
				(*this)[k][rowIndex[j]] = (*this)[k][columnIndex[j]];
				(*this)[k][columnIndex[j]] = d;
			}
		}
	}

	return true;
}

/*
============
CMatXD::inverse_UpdateRankOne

  Updates the in-place inverse using the Sherman-Morrison formula to obtain the inverse for the matrix: A + alpha * v * w'
============
*/
bool CMatXD::inverse_UpdateRankOne( const CVecXD &v, const CVecXD &w, float alpha ) {
	int i, j;
	float beta, s;
	CVecXD y, z;

	SMF_ASSERT( numRows == numColumns );
	SMF_ASSERT( v.getSize() >= numColumns );
	SMF_ASSERT( w.getSize() >= numRows );

	y.setData( numRows, VECX_ALLOCAFLOAT( numRows ) );
	z.setData( numRows, VECX_ALLOCAFLOAT( numRows ) );

	multiply( y, v );
	transposeMultiply( z, w );
	beta = 1.0f + ( w * y );

	if ( beta == 0.0f ) {
		return false;
	}

	alpha /= beta;

	for ( i = 0; i < numRows; i++ ) {
		s = y[i] * alpha;
		for ( j = 0; j < numColumns; j++ ) {
			(*this)[i][j] -= s * z[j];
		}
	}
	return true;
}

/*
============
CMatXD::inverse_UpdateRowColumn

  Updates the in-place inverse to obtain the inverse for the matrix:

      [ 0  a  0 ]
  A + [ d  b  e ]
      [ 0  c  0 ]

  where: a = v[0,r-1], b = v[r], c = v[r+1,numRows-1], d = w[0,r-1], w[r] = 0.0f, e = w[r+1,numColumns-1]
============
*/
bool CMatXD::inverse_UpdateRowColumn( const CVecXD &v, const CVecXD &w, int r ) {
	CVecXD s;

	SMF_ASSERT( numRows == numColumns );
	SMF_ASSERT( v.getSize() >= numColumns );
	SMF_ASSERT( w.getSize() >= numRows );
	SMF_ASSERT( r >= 0 && r < numRows && r < numColumns );
	SMF_ASSERT( w[r] == 0.0f );

	s.setData( MAX( numRows, numColumns ), VECX_ALLOCAFLOAT( MAX( numRows, numColumns ) ) );
	s.toZero();
	s[r] = 1.0f;

	if ( !inverse_UpdateRankOne( v, s, 1.0f ) ) {
		return false;
	}
	if ( !inverse_UpdateRankOne( s, w, 1.0f ) ) {
		return false;
	}
	return true;
}

/*
============
CMatXD::inverse_UpdateIncrement

  Updates the in-place inverse to obtain the inverse for the matrix:

  [ A  a ]
  [ c  b ]

  where: a = v[0,numRows-1], b = v[numRows], c = w[0,numColumns-1], w[numColumns] = 0
============
*/
bool CMatXD::inverse_UpdateIncrement( const CVecXD &v, const CVecXD &w ) {
	CVecXD v2;

	SMF_ASSERT( numRows == numColumns );
	SMF_ASSERT( v.getSize() >= numRows+1 );
	SMF_ASSERT( w.getSize() >= numColumns+1 );

	changeSize( numRows+1, numColumns+1, true );
	(*this)[numRows-1][numRows-1] = 1.0f;

	v2.setData( numRows, VECX_ALLOCAFLOAT( numRows ) );
	v2 = v;
	v2[numRows-1] -= 1.0f;

	return inverse_UpdateRowColumn( v2, w, numRows-1 );
}

/*
============
CMatXD::inverse_UpdateDecrement

  Updates the in-place inverse to obtain the inverse of the matrix with row r and column r removed.
  v and w should store the column and row of the original matrix respectively.
============
*/
bool CMatXD::inverse_UpdateDecrement( const CVecXD &v, const CVecXD &w, int r ) {
	CVecXD v1, w1;

	SMF_ASSERT( numRows == numColumns );
	SMF_ASSERT( v.getSize() >= numRows );
	SMF_ASSERT( w.getSize() >= numColumns );
	SMF_ASSERT( r >= 0 && r < numRows && r < numColumns );

	v1.setData( numRows, VECX_ALLOCAFLOAT( numRows ) );
	w1.setData( numRows, VECX_ALLOCAFLOAT( numRows ) );

	// update the row and column to identity
	v1 = -v;
	w1 = -w;
	v1[r] += 1.0f;
	w1[r] = 0.0f;

	if ( !inverse_UpdateRowColumn( v1, w1, r ) ) {
		return false;
	}

	// physically remove the row and column
	update_Decrement( r );

	return true;
}

/*
============
CMatXD::inverse_solve

  solve Ax = b with A inverted
============
*/
void CMatXD::inverse_solve( CVecXD &x, const CVecXD &b ) const {
	multiply( x, b );
}

/*
============
CMatXD::lu_Factor

  in-place factorization: LU
  L is a triangular matrix stored in the lower triangle.
  L has ones on the diagonal that are not stored.
  U is a triangular matrix stored in the upper triangle.
  If index != NULL partial pivoting is used for numerical stability.
  If index != NULL it must point to an array of numRow integers and is used to keep track of the row permutation.
  If det != NULL the determinant of the matrix is calculated and stored.
============
*/
bool CMatXD::lu_Factor( int *index, float *det ) {
	int i, j, k, newi, min;
	double s, t, d, w;

	// if partial pivoting should be used
	if ( index ) {
		for ( i = 0; i < numRows; i++ ) {
			index[i] = i;
		}
	}

	w = 1.0f;
	min = MIN( numRows, numColumns );
	for ( i = 0; i < min; i++ ) {

		newi = i;
		s = CMath::fabs( (*this)[i][i] );

		if ( index ) {
			// find the largest absolute pivot
			for ( j = i + 1; j < numRows; j++ ) {
				t = CMath::fabs( (*this)[j][i] );
				if ( t > s ) {
					newi = j;
					s = t;
				}
			}
		}

		if ( s == 0.0f ) {
			return false;
		}

		if ( newi != i ) {

			w = -w;

			// swap index elements
			k = index[i];
			index[i] = index[newi];
			index[newi] = k;

			// swap rows
			for ( j = 0; j < numColumns; j++ ) {
				t = (*this)[newi][j];
				(*this)[newi][j] = (*this)[i][j];
				(*this)[i][j] = t;
			}
		}

		if ( i < numRows ) {
			d = 1.0f / (*this)[i][i];
			for ( j = i + 1; j < numRows; j++ ) {
				(*this)[j][i] *= d;
			}
		}

		if ( i < min-1 ) {
			for ( j = i + 1; j < numRows; j++ ) {
				d = (*this)[j][i];
				for ( k = i + 1; k < numColumns; k++ ) {
					(*this)[j][k] -= d * (*this)[i][k];
				}
			}
		}
	}

	if ( det ) {
		for ( i = 0; i < numRows; i++ ) {
			w *= (*this)[i][i];
		}
		*det = w;
	}

	return true;
}

/*
============
CMatXD::lu_UpdateRankOne

  Updates the in-place LU factorization to obtain the factors for the matrix: LU + alpha * v * w'
============
*/
bool CMatXD::lu_UpdateRankOne( const CVecXD &v, const CVecXD &w, float alpha, int *index ) {
	int i, j, max;
	float *y, *z;
	double diag, beta, p0, p1, d;

	SMF_ASSERT( v.getSize() >= numColumns );
	SMF_ASSERT( w.getSize() >= numRows );

	y = (float *) _allocafloat16( v.getSize() * sizeof( float ) );
	z = (float *) _allocafloat16( w.getSize() * sizeof( float ) );

	if ( index != NULL ) {
		for ( i = 0; i < numRows; i++ ) {
			y[i] = alpha * v[index[i]];
		}
	} else {
		for ( i = 0; i < numRows; i++ ) {
			y[i] = alpha * v[i];
		}
	}

	memcpy( z, w.toFloatPtr(), w.getSize() * sizeof( float ) );

	max = MIN( numRows, numColumns );
	for ( i = 0; i < max; i++ ) {
		diag = (*this)[i][i];

		p0 = y[i];
		p1 = z[i];
		diag += p0 * p1;

		if ( diag == 0.0f ) {
			return false;
		}

		beta = p1 / diag;

		(*this)[i][i] = diag;

		for ( j = i+1; j < numColumns; j++ ) {

			d = (*this)[i][j];

			d += p0 * z[j];
			z[j] -= beta * d;

			(*this)[i][j] = d;
		}

		for ( j = i+1; j < numRows; j++ ) {

			d = (*this)[j][i];

			y[j] -= p0 * d;
			d += beta * y[j];

			(*this)[j][i] = d;
		}
	}
	return true;
}

/*
============
CMatXD::lu_UpdateRowColumn

  Updates the in-place LU factorization to obtain the factors for the matrix:

       [ 0  a  0 ]
  LU + [ d  b  e ]
       [ 0  c  0 ]

  where: a = v[0,r-1], b = v[r], c = v[r+1,numRows-1], d = w[0,r-1], w[r] = 0.0f, e = w[r+1,numColumns-1]
============
*/
bool CMatXD::lu_UpdateRowColumn( const CVecXD &v, const CVecXD &w, int r, int *index ) {
#if 0

	CVecXD s;

	SMF_ASSERT( v.getSize() >= numColumns );
	SMF_ASSERT( w.getSize() >= numRows );
	SMF_ASSERT( r >= 0 && r < numRows && r < numColumns );
	SMF_ASSERT( w[r] == 0.0f );

	s.setData( MAX( numRows, numColumns ), VECX_ALLOCA( MAX( numRows, numColumns ) ) );
	s.toZero();
	s[r] = 1.0f;

	if ( !lu_UpdateRankOne( v, s, 1.0f, index ) ) {
		return false;
	}
	if ( !lu_UpdateRankOne( s, w, 1.0f, index ) ) {
		return false;
	}
	return true;

#else

	int i, j, min, max, rp;
	float *y0, *y1, *z0, *z1;
	double diag, beta0, beta1, p0, p1, q0, q1, d;

	SMF_ASSERT( v.getSize() >= numColumns );
	SMF_ASSERT( w.getSize() >= numRows );
	SMF_ASSERT( r >= 0 && r < numColumns && r < numRows );
	SMF_ASSERT( w[r] == 0.0f );

	y0 = (float *) _allocafloat16( v.getSize() * sizeof( float ) );
	z0 = (float *) _allocafloat16( w.getSize() * sizeof( float ) );
	y1 = (float *) _allocafloat16( v.getSize() * sizeof( float ) );
	z1 = (float *) _allocafloat16( w.getSize() * sizeof( float ) );

	if ( index != NULL ) {
		for ( i = 0; i < numRows; i++ ) {
			y0[i] = v[index[i]];
		}
		rp = r;
		for ( i = 0; i < numRows; i++ ) {
			if ( index[i] == r ) {
				rp = i;
				break;
			}
		}
	} else {
		memcpy( y0, v.toFloatPtr(), v.getSize() * sizeof( float ) );
		rp = r;
	}

	memset( y1, 0, v.getSize() * sizeof( float ) );
	y1[rp] = 1.0f;

	memset( z0, 0, w.getSize() * sizeof( float ) );
	z0[r] = 1.0f;

	memcpy( z1, w.toFloatPtr(), w.getSize() * sizeof( float ) );

	// update the beginning of the to be updated row and column
	min = MIN( r, rp );
	for ( i = 0; i < min; i++ ) {
		p0 = y0[i];
		beta1 = z1[i] / (*this)[i][i];

		(*this)[i][r] += p0;
		for ( j = i+1; j < numColumns; j++ ) {
			z1[j] -= beta1 * (*this)[i][j];
		}
		for ( j = i+1; j < numRows; j++ ) {
			y0[j] -= p0 * (*this)[j][i];
		}
		(*this)[rp][i] += beta1;
	}

	// update the lower right corner starting at r,r
	max = MIN( numRows, numColumns );
	for ( i = min; i < max; i++ ) {
		diag = (*this)[i][i];

		p0 = y0[i];
		p1 = z0[i];
		diag += p0 * p1;

		if ( diag == 0.0f ) {
			return false;
		}

		beta0 = p1 / diag;

		q0 = y1[i];
		q1 = z1[i];
		diag += q0 * q1;

		if ( diag == 0.0f ) {
			return false;
		}

		beta1 = q1 / diag;

		(*this)[i][i] = diag;

		for ( j = i+1; j < numColumns; j++ ) {

			d = (*this)[i][j];

			d += p0 * z0[j];
			z0[j] -= beta0 * d;

			d += q0 * z1[j];
			z1[j] -= beta1 * d;

			(*this)[i][j] = d;
		}

		for ( j = i+1; j < numRows; j++ ) {

			d = (*this)[j][i];

			y0[j] -= p0 * d;
			d += beta0 * y0[j];

			y1[j] -= q0 * d;
			d += beta1 * y1[j];

			(*this)[j][i] = d;
		}
	}
	return true;

#endif
}

/*
============
CMatXD::lu_UpdateIncrement

  Updates the in-place LU factorization to obtain the factors for the matrix:

  [ A  a ]
  [ c  b ]

  where: a = v[0,numRows-1], b = v[numRows], c = w[0,numColumns-1], w[numColumns] = 0
============
*/
bool CMatXD::lu_UpdateIncrement( const CVecXD &v, const CVecXD &w, int *index ) {
	int i, j;
	float sum;

	SMF_ASSERT( numRows == numColumns );
	SMF_ASSERT( v.getSize() >= numRows+1 );
	SMF_ASSERT( w.getSize() >= numColumns+1 );

	changeSize( numRows+1, numColumns+1, true );

	// add row to L
	for ( i = 0; i < numRows - 1; i++ ) {
		sum = w[i];
		for ( j = 0; j < i; j++ ) {
			sum -= (*this)[numRows - 1][j] * (*this)[j][i];
		}
		(*this)[numRows - 1 ][i] = sum / (*this)[i][i];
	}

	// add row to the permutation index
	if ( index != NULL ) {
		index[numRows - 1] = numRows - 1;
	}

	// add column to U
	for ( i = 0; i < numRows; i++ ) {
		if ( index != NULL ) {
			sum = v[index[i]];
		} else {
			sum = v[i];
		}
		for ( j = 0; j < i; j++ ) {
			sum -= (*this)[i][j] * (*this)[j][numRows - 1];
		}
		(*this)[i][numRows - 1] = sum;
	}

	return true;
}

/*
============
CMatXD::lu_UpdateDecrement

  Updates the in-place LU factorization to obtain the factors for the matrix with row r and column r removed.
  v and w should store the column and row of the original matrix respectively.
  If index != NULL then u should store row index[r] of the original matrix. If index == NULL then u = w.
============
*/
bool CMatXD::lu_UpdateDecrement( const CVecXD &v, const CVecXD &w, const CVecXD &u, int r, int *index ) {
	int i, p;
	CVecXD v1, w1;

	SMF_ASSERT( numRows == numColumns );
	SMF_ASSERT( v.getSize() >= numColumns );
	SMF_ASSERT( w.getSize() >= numRows );
	SMF_ASSERT( r >= 0 && r < numRows && r < numColumns );

	v1.setData( numRows, VECX_ALLOCAFLOAT( numRows ) );
	w1.setData( numRows, VECX_ALLOCAFLOAT( numRows ) );

	if ( index != NULL ) {

		// find the pivot row
		for ( p = i = 0; i < numRows; i++ ) {
			if ( index[i] == r ) {
				p = i;
				break;
			}
		}

		// update the row and column to identity
		v1 = -v;
		w1 = -u;

		if ( p != r ) {
			swapElements( v1[index[r]], v1[index[p]] );
			swapElements( index[r], index[p] );
		}

		v1[r] += 1.0f;
		w1[r] = 0.0f;

		if ( !lu_UpdateRowColumn( v1, w1, r, index ) ) {
			return false;
		}

		if ( p != r ) {

			if ( CMath::fabs( u[p] ) < 1e-4f ) {
				// NOTE: an additional row interchange is required for numerical stability
			}

			// move row index[r] of the original matrix to row index[p] of the original matrix
			v1.toZero();
			v1[index[p]] = 1.0f;
			w1 = u - w;

			if ( !lu_UpdateRankOne( v1, w1, 1.0f, index ) ) {
				return false;
			}
		}

		// remove the row from the permutation index
		for ( i = r; i < numRows - 1; i++ ) {
			index[i] = index[i+1];
		}
		for ( i = 0; i < numRows - 1; i++ ) {
			if ( index[i] > r ) {
				index[i]--;
			}
		}

	} else {

		v1 = -v;
		w1 = -w;
		v1[r] += 1.0f;
		w1[r] = 0.0f;

		if ( !lu_UpdateRowColumn( v1, w1, r, index ) ) {
			return false;
		}
	}

	// physically remove the row and column
	update_Decrement( r );

	return true;
}

/*
============
CMatXD::lu_solve

  solve Ax = b with A factored in-place as: LU
============
*/
void CMatXD::lu_solve( CVecXD &x, const CVecXD &b, const int *index ) const {
	int i, j;
	double sum;

	SMF_ASSERT( x.getSize() == numColumns && b.getSize() == numRows );

	// solve L
	for ( i = 0; i < numRows; i++ ) {
		if ( index != NULL ) {
			sum = b[index[i]];
		} else {
			sum = b[i];
		}
		for ( j = 0; j < i; j++ ) {
			sum -= (*this)[i][j] * x[j];
		}
		x[i] = sum;
	}

	// solve U
	for ( i = numRows - 1; i >= 0; i-- ) {
		sum = x[i];
		for ( j = i + 1; j < numRows; j++ ) {
			sum -= (*this)[i][j] * x[j];
		}
		x[i] = sum / (*this)[i][i];
	}
}

/*
============
CMatXD::lu_Inverse

  Calculates the inverse of the matrix which is factored in-place as LU
============
*/
void CMatXD::lu_Inverse( CMatXD &inv, const int *index ) const {
	int i, j;
	CVecXD x, b;

	SMF_ASSERT( numRows == numColumns );

	x.setData( numRows, VECX_ALLOCAFLOAT( numRows ) );
	b.setData( numRows, VECX_ALLOCAFLOAT( numRows ) );
	b.toZero();
	inv.setSize( numRows, numColumns );

	for ( i = 0; i < numRows; i++ ) {

		b[i] = 1.0f;
		lu_solve( x, b, index );
		for ( j = 0; j < numRows; j++ ) {
			inv[j][i] = x[j];
		}
		b[i] = 0.0f;
	}
}

/*
============
CMatXD::lu_UnpackFactors

  Unpacks the in-place LU factorization.
============
*/
void CMatXD::lu_UnpackFactors( CMatXD &L, CMatXD &U ) const {
	int i, j;

	L.zero( numRows, numColumns );
	U.zero( numRows, numColumns );
	for ( i = 0; i < numRows; i++ ) {
		for ( j = 0; j < i; j++ ) {
			L[i][j] = (*this)[i][j];
		}
		L[i][i] = 1.0f;
		for ( j = i; j < numColumns; j++ ) {
			U[i][j] = (*this)[i][j];
		}
	}
}

/*
============
CMatXD::lu_MultiplyFactors

  Multiplies the factors of the in-place LU factorization to form the original matrix.
============
*/
void CMatXD::lu_MultiplyFactors( CMatXD &m, const int *index ) const {
	int r, rp, i, j;
	double sum;

	m.setSize( numRows, numColumns );

	for ( r = 0; r < numRows; r++ ) {

		if ( index != NULL ) {
			rp = index[r];
		} else {
			rp = r;
		}

		// calculate row of matrix
		for ( i = 0; i < numColumns; i++ ) {
			if ( i >= r ) {
				sum = (*this)[r][i];
			} else {
				sum = 0.0f;
			}
			for ( j = 0; j <= i && j < r; j++ ) {
				sum += (*this)[r][j] * (*this)[j][i];
			}
			m[rp][i] = sum;
		}
	}
}

/*
============
CMatXD::qr_Factor

  in-place factorization: QR
  Q is an orthogonal matrix represented as a product of Householder matrices stored in the lower triangle and c.
  R is a triangular matrix stored in the upper triangle except for the diagonal elements which are stored in d.
  The initial matrix has to be square.
============
*/
bool CMatXD::qr_Factor( CVecXD &c, CVecXD &d ) {
	int i, j, k;
	double scale, s, t, sum;
	bool singular = false;

	SMF_ASSERT( numRows == numColumns );
	SMF_ASSERT( c.getSize() >= numRows && d.getSize() >= numRows );

	for ( k = 0; k < numRows-1; k++ ) {

		scale = 0.0f;
		for ( i = k; i < numRows; i++ ) {
			s = CMath::fabs( (*this)[i][k] );
			if ( s > scale ) {
				scale = s;
			}
		}
		if ( scale == 0.0f ) {
			singular = true;
			c[k] = d[k] = 0.0f;
		} else {

			s = 1.0f / scale;
			for ( i = k; i < numRows; i++ ) {
				(*this)[i][k] *= s;
			}

			sum = 0.0f;
			for ( i = k; i < numRows; i++ ) {
				s = (*this)[i][k];
				sum += s * s;
			}

			s = CMath::sqrt( sum );
			if ( (*this)[k][k] < 0.0f ) {
				s = -s;
			}
			(*this)[k][k] += s;
			c[k] = s * (*this)[k][k];
			d[k] = -scale * s;

			for ( j = k+1; j < numRows; j++ ) {

				sum = 0.0f;
				for ( i = k; i < numRows; i++ ) {
					sum += (*this)[i][k] * (*this)[i][j];
				}
				t = sum / c[k];
				for ( i = k; i < numRows; i++ ) {
					(*this)[i][j] -= t * (*this)[i][k];
				}
			}
		}
	}
	d[numRows-1] = (*this)[ (numRows-1) ][ (numRows-1) ];
	if ( d[numRows-1] == 0.0f ) {
		singular = true;
	}

	return !singular;
}

/*
============
CMatXD::qr_Rotate

  Performs a Jacobi rotation on the rows i and i+1 of the unpacked QR factors.
============
*/
void CMatXD::qr_Rotate( CMatXD &R, int i, float a, float b ) {
	int j;
	float f, c, s, w, y;

	if ( a == 0.0f ) {
		c = 0.0f;
		s = ( b >= 0.0f ) ? 1.0f : -1.0f;
	} else if ( CMath::fabs( a ) > CMath::fabs( b ) ) {
		f = b / a;
		c = CMath::fabs( 1.0f / CMath::sqrt( 1.0f + f * f ) );
		if ( a < 0.0f ) {
			c = -c;
		}
		s = f * c;
	} else {
		f = a / b;
		s = CMath::fabs( 1.0f / CMath::sqrt( 1.0f + f * f ) );
		if ( b < 0.0f ) {
			s = -s;
		}
		c = f * s;
	}
	for ( j = i; j < numRows; j++ ) {
		y = R[i][j];
		w = R[i+1][j];
		R[i][j] = c * y - s * w;
		R[i+1][j] = s * y + c * w;
	}
	for ( j = 0; j < numRows; j++ ) {
		y = (*this)[j][i];
		w = (*this)[j][i+1];
		(*this)[j][i] = c * y - s * w;
		(*this)[j][i+1] = s * y + c * w;
	}
}

/*
============
CMatXD::qr_UpdateRankOne

  Updates the unpacked QR factorization to obtain the factors for the matrix: QR + alpha * v * w'
============
*/
bool CMatXD::qr_UpdateRankOne( CMatXD &R, const CVecXD &v, const CVecXD &w, float alpha ) {
	int i, k;
	float f;
	CVecXD u;

	SMF_ASSERT( v.getSize() >= numColumns );
	SMF_ASSERT( w.getSize() >= numRows );

	u.setData( v.getSize(), VECX_ALLOCAFLOAT( v.getSize() ) );
	transposeMultiply( u, v );
	u *= alpha;

	for ( k = v.getSize()-1; k > 0; k-- ) {
		if ( u[k] != 0.0f ) {
			break;
		}
	}
	for ( i = k-1; i >= 0; i-- ) {
		qr_Rotate( R, i, u[i], -u[i+1] );
		if ( u[i] == 0.0f ) {
			u[i] = CMath::fabs( u[i+1] );
		} else if ( CMath::fabs( u[i] ) > CMath::fabs( u[i+1] ) ) {
			f = u[i+1] / u[i];
			u[i] = CMath::fabs( u[i] ) * CMath::sqrt( 1.0f + f * f );
		} else {
			f = u[i] / u[i+1];
			u[i] = CMath::fabs( u[i+1] ) * CMath::sqrt( 1.0f + f * f );
		}
	}
	for ( i = 0; i < v.getSize(); i++ ) {
		R[0][i] += u[0] * w[i];
	}
	for ( i = 0; i < k; i++ ) {
		qr_Rotate( R, i, -R[i][i], R[i+1][i] );
	}
	return true;
}

/*
============
CMatXD::qr_UpdateRowColumn

  Updates the unpacked QR factorization to obtain the factors for the matrix:

       [ 0  a  0 ]
  QR + [ d  b  e ]
       [ 0  c  0 ]

  where: a = v[0,r-1], b = v[r], c = v[r+1,numRows-1], d = w[0,r-1], w[r] = 0.0f, e = w[r+1,numColumns-1]
============
*/
bool CMatXD::qr_UpdateRowColumn( CMatXD &R, const CVecXD &v, const CVecXD &w, int r ) {
	CVecXD s;

	SMF_ASSERT( v.getSize() >= numColumns );
	SMF_ASSERT( w.getSize() >= numRows );
	SMF_ASSERT( r >= 0 && r < numRows && r < numColumns );
	SMF_ASSERT( w[r] == 0.0f );

	s.setData( MAX( numRows, numColumns ), VECX_ALLOCAFLOAT( MAX( numRows, numColumns ) ) );
	s.toZero();
	s[r] = 1.0f;

	if ( !qr_UpdateRankOne( R, v, s, 1.0f ) ) {
		return false;
	}
	if ( !qr_UpdateRankOne( R, s, w, 1.0f ) ) {
		return false;
	}
	return true;
}

/*
============
CMatXD::qr_UpdateIncrement

  Updates the unpacked QR factorization to obtain the factors for the matrix:

  [ A  a ]
  [ c  b ]

  where: a = v[0,numRows-1], b = v[numRows], c = w[0,numColumns-1], w[numColumns] = 0
============
*/
bool CMatXD::qr_UpdateIncrement( CMatXD &R, const CVecXD &v, const CVecXD &w ) {
	CVecXD v2;

	SMF_ASSERT( numRows == numColumns );
	SMF_ASSERT( v.getSize() >= numRows+1 );
	SMF_ASSERT( w.getSize() >= numColumns+1 );

	changeSize( numRows+1, numColumns+1, true );
	(*this)[numRows-1][numRows-1] = 1.0f;

	R.changeSize( R.numRows+1, R.numColumns+1, true );
	R[R.numRows-1][R.numRows-1] = 1.0f;

	v2.setData( numRows, VECX_ALLOCAFLOAT( numRows ) );
	v2 = v;
	v2[numRows-1] -= 1.0f;

	return qr_UpdateRowColumn( R, v2, w, numRows-1 );
}

/*
============
CMatXD::qr_UpdateDecrement

  Updates the unpacked QR factorization to obtain the factors for the matrix with row r and column r removed.
  v and w should store the column and row of the original matrix respectively.
============
*/
bool CMatXD::qr_UpdateDecrement( CMatXD &R, const CVecXD &v, const CVecXD &w, int r ) {
	CVecXD v1, w1;

	SMF_ASSERT( numRows == numColumns );
	SMF_ASSERT( v.getSize() >= numRows );
	SMF_ASSERT( w.getSize() >= numColumns );
	SMF_ASSERT( r >= 0 && r < numRows && r < numColumns );

	v1.setData( numRows, VECX_ALLOCAFLOAT( numRows ) );
	w1.setData( numRows, VECX_ALLOCAFLOAT( numRows ) );

	// update the row and column to identity
	v1 = -v;
	w1 = -w;
	v1[r] += 1.0f;
	w1[r] = 0.0f;

	if ( !qr_UpdateRowColumn( R, v1, w1, r ) ) {
		return false;
	}

	// physically remove the row and column
	update_Decrement( r );
	R.update_Decrement( r );

	return true;
}

/*
============
CMatXD::qr_solve

  solve Ax = b with A factored in-place as: QR
============
*/
void CMatXD::qr_solve( CVecXD &x, const CVecXD &b, const CVecXD &c, const CVecXD &d ) const {
	int i, j;
	double sum, t;

	SMF_ASSERT( numRows == numColumns );
	SMF_ASSERT( x.getSize() >= numRows && b.getSize() >= numRows );
	SMF_ASSERT( c.getSize() >= numRows && d.getSize() >= numRows );

	for ( i = 0; i < numRows; i++ ) {
		x[i] = b[i];
	}

	// multiply b with transpose of Q
	for ( i = 0; i < numRows-1; i++ ) {

		sum = 0.0f;
		for ( j = i; j < numRows; j++ ) {
			sum += (*this)[j][i] * x[j];
		}
		t = sum / c[i];
		for ( j = i; j < numRows; j++ ) {
			x[j] -= t * (*this)[j][i];
		}
	}

	// backsubstitution with R
	for ( i = numRows-1; i >= 0; i-- ) {

		sum = x[i];
		for ( j = i + 1; j < numRows; j++ ) {
			sum -= (*this)[i][j] * x[j];
		}
		x[i] = sum / d[i];
	}
}

/*
============
CMatXD::qr_solve

  solve Ax = b with A factored as: QR
============
*/
void CMatXD::qr_solve( CVecXD &x, const CVecXD &b, const CMatXD &R ) const {
	int i, j;
	double sum;

	SMF_ASSERT( numRows == numColumns );

	// multiply b with transpose of Q
	transposeMultiply( x, b );

	// backsubstitution with R
	for ( i = numRows-1; i >= 0; i-- ) {

		sum = x[i];
		for ( j = i + 1; j < numRows; j++ ) {
			sum -= R[i][j] * x[j];
		}
		x[i] = sum / R[i][i];
	}
}

/*
============
CMatXD::qr_Inverse

  Calculates the inverse of the matrix which is factored in-place as: QR
============
*/
void CMatXD::qr_Inverse( CMatXD &inv, const CVecXD &c, const CVecXD &d ) const {
	int i, j;
	CVecXD x, b;

	SMF_ASSERT( numRows == numColumns );

	x.setData( numRows, VECX_ALLOCAFLOAT( numRows ) );
	b.setData( numRows, VECX_ALLOCAFLOAT( numRows ) );
	b.toZero();
	inv.setSize( numRows, numColumns );

	for ( i = 0; i < numRows; i++ ) {

		b[i] = 1.0f;
		qr_solve( x, b, c, d );
		for ( j = 0; j < numRows; j++ ) {
			inv[j][i] = x[j];
		}
		b[i] = 0.0f;
	}
}

/*
============
CMatXD::qr_UnpackFactors

  Unpacks the in-place QR factorization.
============
*/
void CMatXD::qr_UnpackFactors( CMatXD &Q, CMatXD &R, const CVecXD &c, const CVecXD &d ) const {
	int i, j, k;
	double sum;

	Q.identity( numRows, numColumns );
	for ( i = 0; i < numColumns-1; i++ ) {
		if ( c[i] == 0.0f ) {
			continue;
		}
		for ( j = 0; j < numRows; j++ ) {
			sum = 0.0f;
			for ( k = i; k < numColumns; k++ ) {
				sum += (*this)[k][i] * Q[j][k];
			}
			sum /= c[i];
			for ( k = i; k < numColumns; k++ ) {
				Q[j][k] -= sum * (*this)[k][i];
			}
		}
	}

	R.zero( numRows, numColumns );
	for ( i = 0; i < numRows; i++ ) {
		R[i][i] = d[i];
		for ( j = i+1; j < numColumns; j++ ) {
			R[i][j] = (*this)[i][j];
		}
	}
}

/*
============
CMatXD::qr_MultiplyFactors

  Multiplies the factors of the in-place QR factorization to form the original matrix.
============
*/
void CMatXD::qr_MultiplyFactors( CMatXD &m, const CVecXD &c, const CVecXD &d ) const {
	int i, j, k;
	double sum;
	CMatXD Q;

	Q.identity( numRows, numColumns );
	for ( i = 0; i < numColumns-1; i++ ) {
		if ( c[i] == 0.0f ) {
			continue;
		}
		for ( j = 0; j < numRows; j++ ) {
			sum = 0.0f;
			for ( k = i; k < numColumns; k++ ) {
				sum += (*this)[k][i] * Q[j][k];
			}
			sum /= c[i];
			for ( k = i; k < numColumns; k++ ) {
				Q[j][k] -= sum * (*this)[k][i];
			}
		}
	}

	for ( i = 0; i < numRows; i++ ) {
		for ( j = 0; j < numColumns; j++ ) {
			sum = Q[i][j] * d[i];
			for ( k = 0; k < i; k++ ) {
				sum += Q[i][k] * (*this)[j][k];
			}
			m[i][j] = sum;
		}
	}
}

/*
============
CMatXD::Pythag

  Computes (a^2 + b^2)^1/2 without underflow or overflow.
============
*/
float CMatXD::Pythag( float a, float b ) const {
	double at, bt, ct;

	at = CMath::fabs( a );
	bt = CMath::fabs( b );
	if ( at > bt ) {
		ct = bt / at;
		return at * CMath::sqrt( 1.0f + ct * ct );
	} else {
		if ( bt ) {
			ct = at / bt;
			return bt * CMath::sqrt( 1.0f + ct * ct );
		} else {
			return 0.0f;
		}
	}
}

/*
============
CMatXD::svd_BiDiag
============
*/
void CMatXD::svd_BiDiag( CVecXD &w, CVecXD &rv1, float &anorm ) {
	int i, j, k, l;
	double f, h, r, g, s, scale;

	anorm = 0.0f;
	g = s = scale = 0.0f;
	for ( i = 0; i < numColumns; i++ ) {
		l = i + 1;
		rv1[i] = scale * g;
		g = s = scale = 0.0f;
		if ( i < numRows ) {
			for ( k = i; k < numRows; k++ ) {
				scale += CMath::fabs( (*this)[k][i] );
			}
			if ( scale ) {
				for ( k = i; k < numRows; k++ ) {
					(*this)[k][i] /= scale;
					s += (*this)[k][i] * (*this)[k][i];
				}
				f = (*this)[i][i];
				g = CMath::sqrt( s );
				if ( f >= 0.0f ) {
					g = -g;
				}
				h = f * g - s;
				(*this)[i][i] = f - g;
				if ( i != (numColumns-1) ) {
					for ( j = l; j < numColumns; j++ ) {
						for ( s = 0.0f, k = i; k < numRows; k++ ) {
							s += (*this)[k][i] * (*this)[k][j];
						}
						f = s / h;
						for ( k = i; k < numRows; k++ ) {
							(*this)[k][j] += f * (*this)[k][i];
						}
					}
				}
				for ( k = i; k < numRows; k++ ) {
					(*this)[k][i] *= scale;
				}
			}
		}
		w[i] = scale * g;
		g = s = scale = 0.0f;
		if ( i < numRows && i != (numColumns-1) ) {
			for ( k = l; k < numColumns; k++ ) {
				scale += CMath::fabs( (*this)[i][k] );
			}
			if ( scale ) {
				for ( k = l; k < numColumns; k++ ) {
					(*this)[i][k] /= scale;
					s += (*this)[i][k] * (*this)[i][k];
				}
				f = (*this)[i][l];
				g = CMath::sqrt( s );
				if ( f >= 0.0f ) {
					g = -g;
				}
				h = 1.0f / ( f * g - s );
				(*this)[i][l] = f - g;
				for ( k = l; k < numColumns; k++ ) {
					rv1[k] = (*this)[i][k] * h;
				}
				if ( i != (numRows-1) ) {
					for ( j = l; j < numRows; j++ ) {
						for ( s = 0.0f, k = l; k < numColumns; k++ ) {
							s += (*this)[j][k] * (*this)[i][k];
						}
						for ( k = l; k < numColumns; k++ ) {
							(*this)[j][k] += s * rv1[k];
						}
					}
				}
				for ( k = l; k < numColumns; k++ ) {
					(*this)[i][k] *= scale;
				}
			}
		}
		r = CMath::fabs( w[i] ) + CMath::fabs( rv1[i] );
		if ( r > anorm ) {
			anorm = r;
		}
	}
}

/*
============
CMatXD::svd_InitialWV
============
*/
void CMatXD::svd_InitialWV( CVecXD &w, CMatXD &V, CVecXD &rv1 ) {
	int i, j, k, l;
	double f, g, s;

	g = 0.0f;
	for ( i = (numColumns-1); i >= 0; i-- ) {
		l = i + 1;
		if ( i < ( numColumns - 1 ) ) {
			if ( g ) {
				for ( j = l; j < numColumns; j++ ) {
					V[j][i] = ((*this)[i][j] / (*this)[i][l]) / g;
				}
				// double division to reduce underflow
				for ( j = l; j < numColumns; j++ ) {
					for ( s = 0.0f, k = l; k < numColumns; k++ ) {
						s += (*this)[i][k] * V[k][j];
					}
					for ( k = l; k < numColumns; k++ ) {
						V[k][j] += s * V[k][i];
					}
				}
			}
			for ( j = l; j < numColumns; j++ ) {
				V[i][j] = V[j][i] = 0.0f;
			}
		}
		V[i][i] = 1.0f;
		g = rv1[i];
	}
	for ( i = numColumns - 1 ; i >= 0; i-- ) {
		l = i + 1;
		g = w[i];
		if ( i < (numColumns-1) ) {
			for ( j = l; j < numColumns; j++ ) {
				(*this)[i][j] = 0.0f;
			}
		}
		if ( g ) {
			g = 1.0f / g;
			if ( i != (numColumns-1) ) {
				for ( j = l; j < numColumns; j++ ) {
					for ( s = 0.0f, k = l; k < numRows; k++ ) {
						s += (*this)[k][i] * (*this)[k][j];
					}
					f = (s / (*this)[i][i]) * g;
					for ( k = i; k < numRows; k++ ) {
						(*this)[k][j] += f * (*this)[k][i];
					}
				}
			}
			for ( j = i; j < numRows; j++ ) {
				(*this)[j][i] *= g;
			}
		}
		else {
			for ( j = i; j < numRows; j++ ) {
				(*this)[j][i] = 0.0f;
			}
		}
		(*this)[i][i] += 1.0f;
	}
}

/*
============
CMatXD::svd_Factor

  in-place factorization: U * diag(w) * V.transpose()
  known as the Singular Value Decomposition.
  U is a column-orthogonal matrix which overwrites the original matrix.
  w is a diagonal matrix with all elements >= 0 which are the singular values.
  V is the transpose of an orthogonal matrix.
============
*/
bool CMatXD::svd_Factor( CVecXD &w, CMatXD &V ) {
	int flag, i, its, j, jj, k, l, nm;
	double c, f, h, s, x, y, z, r, g = 0.0f;
	float anorm = 0.0f;
	CVecXD rv1;

	if ( numRows < numColumns ) {
		return false;
	}

	rv1.setData( numColumns, VECX_ALLOCAFLOAT( numColumns ) );
	rv1.toZero();
	w.zero( numColumns );
	V.zero( numColumns, numColumns );

	svd_BiDiag( w, rv1, anorm );
	svd_InitialWV( w, V, rv1 );

	for ( k = numColumns - 1; k >= 0; k-- ) {
		for ( its = 1; its <= 30; its++ ) {
			flag = 1;
			nm = 0;
			for ( l = k; l >= 0; l-- ) {
				nm = l - 1;
				if ( ( CMath::fabs( rv1[l] ) + anorm ) == anorm /* CMath::fabs( rv1[l] ) < CMath::FLT_EPSILON */ ) {
					flag = 0;
					break;
				}
				if ( ( CMath::fabs( w[nm] ) + anorm ) == anorm /* CMath::fabs( w[nm] ) < CMath::FLT_EPSILON */ ) {
					break;
				}
			}
			if ( flag ) {
				c = 0.0f;
				s = 1.0f;
				for ( i = l; i <= k; i++ ) {
					f = s * rv1[i];

					if ( ( CMath::fabs( f ) + anorm ) != anorm /* CMath::fabs( f ) > CMath::FLT_EPSILON */ ) {
						g = w[i];
						h = Pythag( f, g );
						w[i] = h;
						h = 1.0f / h;
						c = g * h;
						s = -f * h;
						for ( j = 0; j < numRows; j++ ) {
							y = (*this)[j][nm];
							z = (*this)[j][i];
							(*this)[j][nm] = y * c + z * s;
							(*this)[j][i] = z * c - y * s;
						}
					}
				}
			}
			z = w[k];
			if ( l == k ) {
				if ( z < 0.0f ) {
					w[k] = -z;
					for ( j = 0; j < numColumns; j++ ) {
						V[j][k] = -V[j][k];
					}
				}
				break;
			}
			if ( its == 30 ) {
				return false;		// no convergence
			}
			x = w[l];
			nm = k - 1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ( ( y - z ) * ( y + z ) + ( g - h ) * ( g + h ) ) / ( 2.0f * h * y );
			g = Pythag( f, 1.0f );
			r = ( f >= 0.0f ? g : - g );
			f= ( ( x - z ) * ( x + z ) + h * ( ( y / ( f + r ) ) - h ) ) / x;
			c = s = 1.0f;
			for ( j = l; j <= nm; j++ ) {
				i = j + 1;
				g = rv1[i];
				y = w[i];
				h = s * g;
				g = c * g;
				z = Pythag( f, h );
				rv1[j] = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y = y * c;
				for ( jj = 0; jj < numColumns; jj++ ) {
					x = V[jj][j];
					z = V[jj][i];
					V[jj][j] = x * c + z * s;
					V[jj][i] = z * c - x * s;
				}
				z = Pythag( f, h );
				w[j] = z;
				if ( z ) {
					z = 1.0f / z;
					c = f * z;
					s = h * z;
				}
				f = ( c * g ) + ( s * y );
				x = ( c * y ) - ( s * g );
				for ( jj = 0; jj < numRows; jj++ ) {
					y = (*this)[jj][j];
					z = (*this)[jj][i];
					(*this)[jj][j] = y * c + z * s;
					(*this)[jj][i] = z * c - y * s;
				}
			}
			rv1[l] = 0.0f;
			rv1[k] = f;
			w[k] = x;
		}
	}
	return true;
}

/*
============
CMatXD::svd_solve

  solve Ax = b with A factored as: U * diag(w) * V.transpose()
============
*/
void CMatXD::svd_solve( CVecXD &x, const CVecXD &b, const CVecXD &w, const CMatXD &V ) const {
	int i, j;
	double sum;
	CVecXD tmp;

	SMF_ASSERT( x.getSize() >= numColumns );
	SMF_ASSERT( b.getSize() >= numColumns );
	SMF_ASSERT( w.getSize() == numColumns );
	SMF_ASSERT( V.getNumRows() == numColumns && V.getNumColumns() == numColumns );

	tmp.setData( numColumns, VECX_ALLOCAFLOAT( numColumns ) );

	for ( i = 0; i < numColumns; i++ ) {
		sum = 0.0f;
		if ( w[i] >= CMath::FLT_EPSILON ) {
			for ( j = 0; j < numRows; j++ ) {
				sum += (*this)[j][i] * b[j];
			}
			sum /= w[i];
		}
		tmp[i] = sum;
	}
	for ( i = 0; i < numColumns; i++ ) {
		sum = 0.0f;
		for ( j = 0; j < numColumns; j++ ) {
			sum += V[i][j] * tmp[j];
		}
		x[i] = sum;
	}
}

/*
============
CMatXD::svd_Inverse

  Calculates the inverse of the matrix which is factored in-place as: U * diag(w) * V.transpose()
============
*/
void CMatXD::svd_Inverse( CMatXD &inv, const CVecXD &w, const CMatXD &V ) const {
	int i, j, k;
	double wi, sum;
	CMatXD V2;

	SMF_ASSERT( numRows == numColumns );

	V2 = V;

	// V * [diag(1/w[i])]
	for ( i = 0; i < numRows; i++ ) {
		wi = w[i];
		wi = ( wi < CMath::FLT_EPSILON ) ? 0.0f : 1.0f / wi;
		for ( j = 0; j < numColumns; j++ ) {
			V2[j][i] *= wi;
		}
	}

	// V * [diag(1/w[i])] * Ut
	for ( i = 0; i < numRows; i++ ) {
		for ( j = 0; j < numColumns; j++ ) {
			sum = V2[i][0] * (*this)[j][0];
			for ( k = 1; k < numColumns; k++ ) {
				sum += V2[i][k] * (*this)[j][k];
			}
			inv[i][j] = sum;
		}
	}
}

/*
============
CMatXD::svd_MultiplyFactors

  Multiplies the factors of the in-place SVD factorization to form the original matrix.
============
*/
void CMatXD::svd_MultiplyFactors( CMatXD &m, const CVecXD &w, const CMatXD &V ) const {
	int r, i, j;
	double sum;

	m.setSize( numRows, V.getNumRows() );

	for ( r = 0; r < numRows; r++ ) {
		// calculate row of matrix
		if ( w[r] >= CMath::FLT_EPSILON ) {
			for ( i = 0; i < V.getNumRows(); i++ ) {
				sum = 0.0f;
				for ( j = 0; j < numColumns; j++ ) {
					sum += (*this)[r][j] * V[i][j];
				}
				m[r][i] = sum * w[r];
			}
		} else {
			for ( i = 0; i < V.getNumRows(); i++ ) {
				m[r][i] = 0.0f;
			}
		}
	}
}

/*
============
CMatXD::cholesky_Factor

  in-place Cholesky factorization: LL'
  L is a triangular matrix stored in the lower triangle.
  The upper triangle is not cleared.
  The initial matrix has to be symmetric positive definite.
============
*/
bool CMatXD::cholesky_Factor() {
	int i, j, k;
	float *invSqrt;
	double sum;

	SMF_ASSERT( numRows == numColumns );

	invSqrt = (float *) _allocafloat16( numRows * sizeof( float ) );

	for ( i = 0; i < numRows; i++ ) {

		for ( j = 0; j < i; j++ ) {

			sum = (*this)[i][j];
			for ( k = 0; k < j; k++ ) {
				sum -= (*this)[i][k] * (*this)[j][k];
			}
			(*this)[i][j] = sum * invSqrt[j];
		}

		sum = (*this)[i][i];
		for ( k = 0; k < i; k++ ) {
			sum -= (*this)[i][k] * (*this)[i][k];
		}

		if ( sum <= 0.0f ) {
			return false;
		}

		invSqrt[i] = CMath::invSqrt( sum );
		(*this)[i][i] = invSqrt[i] * sum;
	}
	return true;
}

/*
============
CMatXD::cholesky_UpdateRankOne

  Updates the in-place Cholesky factorization to obtain the factors for the matrix: LL' + alpha * v * v'
  If offset > 0 only the lower right corner starting at (offset, offset) is updated.
============
*/
bool CMatXD::cholesky_UpdateRankOne( const CVecXD &v, float alpha, int offset ) {
	int i, j;
	float *y;
	double diag, invDiag, diagSqr, newDiag, newDiagSqr, beta, p, d;

	SMF_ASSERT( numRows == numColumns );
	SMF_ASSERT( v.getSize() >= numRows );
	SMF_ASSERT( offset >= 0 && offset < numRows );

	y = (float *) _allocafloat16( v.getSize() * sizeof( float ) );
	memcpy( y, v.toFloatPtr(), v.getSize() * sizeof( float ) );

	for ( i = offset; i < numColumns; i++ ) {
		p = y[i];
		diag = (*this)[i][i];
		invDiag = 1.0f / diag;
		diagSqr = diag * diag;
		newDiagSqr = diagSqr + alpha * p * p;

		if ( newDiagSqr <= 0.0f ) {
			return false;
		}

		(*this)[i][i] = newDiag = CMath::sqrt( newDiagSqr );

		alpha /= newDiagSqr;
		beta = p * alpha;
		alpha *= diagSqr;

		for ( j = i+1; j < numRows; j++ ) {

			d = (*this)[j][i] * invDiag;

			y[j] -= p * d;
			d += beta * y[j];

			(*this)[j][i] = d * newDiag;
		}
	}
	return true;
}

/*
============
CMatXD::cholesky_UpdateRowColumn

  Updates the in-place Cholesky factorization to obtain the factors for the matrix:

        [ 0  a  0 ]
  LL' + [ a  b  c ]
        [ 0  c  0 ]

  where: a = v[0,r-1], b = v[r], c = v[r+1,numRows-1]
============
*/
bool CMatXD::cholesky_UpdateRowColumn( const CVecXD &v, int r ) {
	int i, j;
	double sum;
	float *original, *y;
	CVecXD addSub;

	SMF_ASSERT( numRows == numColumns );
	SMF_ASSERT( v.getSize() >= numRows );
	SMF_ASSERT( r >= 0 && r < numRows );

	addSub.setData( numColumns, (float *) _allocafloat16( numColumns * sizeof( float ) ) );

	if ( r == 0 ) {

		if ( numColumns == 1 ) {
			double v0 = v[0];
			sum = (*this)[0][0];
			sum = sum * sum;
			sum = sum + v0;
			if ( sum <= 0.0f ) {
				return false;
			}
			(*this)[0][0] = CMath::sqrt( sum );
			return true;
		}
		for ( i = 0; i < numColumns; i++ ) {
			addSub[i] = v[i];
		}

	} else {

		original = (float *) _allocafloat16( numColumns * sizeof( float ) );
		y = (float *) _allocafloat16( numColumns * sizeof( float ) );

		// calculate original row/column of matrix
		for ( i = 0; i < numRows; i++ ) {
			sum = 0.0f;
			for ( j = 0; j <= i; j++ ) {
				sum += (*this)[r][j] * (*this)[i][j];
			}
			original[i] = sum;
		}

		// solve for y in L * y = original + v
		for ( i = 0; i < r; i++ ) {
			sum = original[i] + v[i];
			for ( j = 0; j < i; j++ ) {
				sum -= (*this)[r][j] * (*this)[i][j];
			}
			(*this)[r][i] = sum / (*this)[i][i];
		}

		// if the last row/column of the matrix is updated
		if ( r == numColumns - 1 ) {
			// only calculate new diagonal
			sum = original[r] + v[r];
			for ( j = 0; j < r; j++) {
				sum -= (*this)[r][j] * (*this)[r][j];
			}
			if ( sum <= 0.0f ) {
				return false;
			}
			(*this)[r][r] = CMath::sqrt( sum );
			return true;
		}

		// calculate the row/column to be added to the lower right sub matrix starting at (r, r)
		for ( i = r; i < numColumns; i++ ) {
			sum = 0.0f;
			for ( j = 0; j <= r; j++ ) {
				sum += (*this)[r][j] * (*this)[i][j];
			}
			addSub[i] = v[i] - ( sum - original[i] );
		}
	}

	// add row/column to the lower right sub matrix starting at (r, r)

#if 0

	CVecXD v1, v2;
	double d;

	v1.setData( numColumns, (float *) _alloca16( numColumns * sizeof( float ) ) );
	v2.setData( numColumns, (float *) _alloca16( numColumns * sizeof( float ) ) );

	d = CMath::SQRT_1OVER2;
	v1[r] = ( 0.5f * addSub[r] + 1.0f ) * d;
	v2[r] = ( 0.5f * addSub[r] - 1.0f ) * d;
	for ( i = r+1; i < numColumns; i++ ) {
		v1[i] = v2[i] = addSub[i] * d;
	}

	// update
	if ( !cholesky_UpdateRankOne( v1, 1.0f, r ) ) {
		return false;
	}
	// downdate
	if ( !cholesky_UpdateRankOne( v2, -1.0f, r ) ) {
		return false;
	}

#else

	float *v1, *v2;
	double diag, invDiag, diagSqr, newDiag, newDiagSqr;
	double alpha1, alpha2, beta1, beta2, p1, p2, d;

	v1 = (float *) _allocafloat16( numColumns * sizeof( float ) );
	v2 = (float *) _allocafloat16( numColumns * sizeof( float ) );

	d = CMath::SQRT_1OVER2;
	v1[r] = ( 0.5f * addSub[r] + 1.0f ) * d;
	v2[r] = ( 0.5f * addSub[r] - 1.0f ) * d;
	for ( i = r+1; i < numColumns; i++ ) {
		v1[i] = v2[i] = addSub[i] * d;
	}

	alpha1 = 1.0f;
	alpha2 = -1.0f;

	// simultaneous update/downdate of the sub matrix starting at (r, r)
	for ( i = r; i < numColumns; i++ ) {
		p1 = v1[i];
		diag = (*this)[i][i];
		invDiag = 1.0f / diag;
		diagSqr = diag * diag;
		newDiagSqr = diagSqr + alpha1 * p1 * p1;

		if ( newDiagSqr <= 0.0f ) {
			return false;
		}

		alpha1 /= newDiagSqr;
		beta1 = p1 * alpha1;
		alpha1 *= diagSqr;

		p2 = v2[i];
		diagSqr = newDiagSqr;
		newDiagSqr = diagSqr + alpha2 * p2 * p2;

		if ( newDiagSqr <= 0.0f ) {
			return false;
		}

		(*this)[i][i] = newDiag = CMath::sqrt( newDiagSqr );

		alpha2 /= newDiagSqr;
		beta2 = p2 * alpha2;
		alpha2 *= diagSqr;

		for ( j = i+1; j < numRows; j++ ) {

			d = (*this)[j][i] * invDiag;

			v1[j] -= p1 * d;
			d += beta1 * v1[j];

			v2[j] -= p2 * d;
			d += beta2 * v2[j];

			(*this)[j][i] = d * newDiag;
		}
	}

#endif

	return true;
}

/*
============
CMatXD::cholesky_UpdateIncrement

  Updates the in-place Cholesky factorization to obtain the factors for the matrix:

  [ A  a ]
  [ a  b ]

  where: a = v[0,numRows-1], b = v[numRows]
============
*/
bool CMatXD::cholesky_UpdateIncrement( const CVecXD &v ) {
	int i, j;
	float *x;
	double sum;

	SMF_ASSERT( numRows == numColumns );
	SMF_ASSERT( v.getSize() >= numRows+1 );

	changeSize( numRows+1, numColumns+1, false );

	x = (float *) _allocafloat16( numRows * sizeof( float ) );

	// solve for x in L * x = v
	for ( i = 0; i < numRows - 1; i++ ) {
		sum = v[i];
		for ( j = 0; j < i; j++ ) {
			sum -= (*this)[i][j] * x[j];
		}
		x[i] = sum / (*this)[i][i];
	}

	// calculate new row of L and calculate the square of the diagonal entry
	sum = v[numRows - 1];
	for ( i = 0; i < numRows - 1; i++ ) {
		(*this)[numRows - 1][i] = x[i];
		sum -= x[i] * x[i];
	}

	if ( sum <= 0.0f ) {
		return false;
	}

	// store the diagonal entry
	(*this)[numRows - 1][numRows - 1] = CMath::sqrt( sum );

	return true;
}

/*
============
CMatXD::cholesky_UpdateDecrement

  Updates the in-place Cholesky factorization to obtain the factors for the matrix with row r and column r removed.
  v should store the row of the original matrix.
============
*/
bool CMatXD::cholesky_UpdateDecrement( const CVecXD &v, int r ) {
	CVecXD v1;

	SMF_ASSERT( numRows == numColumns );
	SMF_ASSERT( v.getSize() >= numRows );
	SMF_ASSERT( r >= 0 && r < numRows );

	v1.setData( numRows, VECX_ALLOCAFLOAT( numRows ) );

	// update the row and column to identity
	v1 = -v;
	v1[r] += 1.0f;

	// NOTE:	msvc compiler bug: the this pointer stored in edi is expected to stay
	//			untouched when calling cholesky_UpdateRowColumn in the if statement
#if 0
	if ( !cholesky_UpdateRowColumn( v1, r ) ) {
#else
	bool ret = cholesky_UpdateRowColumn( v1, r );
	if ( !ret ) {
#endif
		return false;
	}

	// physically remove the row and column
	update_Decrement( r );

	return true;
}

/*
============
CMatXD::cholesky_solve

  solve Ax = b with A factored in-place as: LL'
============
*/
void CMatXD::cholesky_solve( CVecXD &x, const CVecXD &b ) const {
	int i, j;
	double sum;

	SMF_ASSERT( numRows == numColumns );
	SMF_ASSERT( x.getSize() >= numRows && b.getSize() >= numRows );

	// solve L
	for ( i = 0; i < numRows; i++ ) {
		sum = b[i];
		for ( j = 0; j < i; j++ ) {
			sum -= (*this)[i][j] * x[j];
		}
		x[i] = sum / (*this)[i][i];
	}

	// solve Lt
	for ( i = numRows - 1; i >= 0; i-- ) {
		sum = x[i];
		for ( j = i + 1; j < numRows; j++ ) {
			sum -= (*this)[j][i] * x[j];
		}
		x[i] = sum / (*this)[i][i];
	}
}

/*
============
CMatXD::cholesky_Inverse

  Calculates the inverse of the matrix which is factored in-place as: LL'
============
*/
void CMatXD::cholesky_Inverse( CMatXD &inv ) const {
	int i, j;
	CVecXD x, b;

	SMF_ASSERT( numRows == numColumns );

	x.setData( numRows, VECX_ALLOCAFLOAT( numRows ) );
	b.setData( numRows, VECX_ALLOCAFLOAT( numRows ) );
	b.toZero();
	inv.setSize( numRows, numColumns );

	for ( i = 0; i < numRows; i++ ) {

		b[i] = 1.0f;
		cholesky_solve( x, b );
		for ( j = 0; j < numRows; j++ ) {
			inv[j][i] = x[j];
		}
		b[i] = 0.0f;
	}
}

/*
============
CMatXD::cholesky_MultiplyFactors

  Multiplies the factors of the in-place Cholesky factorization to form the original matrix.
============
*/
void CMatXD::cholesky_MultiplyFactors( CMatXD &m ) const {
	int r, i, j;
	double sum;

	m.setSize( numRows, numColumns );

	for ( r = 0; r < numRows; r++ ) {

		// calculate row of matrix
		for ( i = 0; i < numRows; i++ ) {
			sum = 0.0f;
			for ( j = 0; j <= i && j <= r; j++ ) {
				sum += (*this)[r][j] * (*this)[i][j];
			}
			m[r][i] = sum;
		}
	}
}

/*
============
CMatXD::ldlt_Factor

  in-place factorization: LDL'
  L is a triangular matrix stored in the lower triangle.
  L has ones on the diagonal that are not stored.
  D is a diagonal matrix stored on the diagonal.
  The upper triangle is not cleared.
  The initial matrix has to be symmetric.
============
*/
bool CMatXD::ldlt_Factor() {
	int i, j, k;
	float *v;
	double d, sum;

	SMF_ASSERT( numRows == numColumns );

	v = (float *) _allocafloat16( numRows * sizeof( float ) );

	for ( i = 0; i < numRows; i++ ) {

		sum = (*this)[i][i];
		for ( j = 0; j < i; j++ ) {
			d = (*this)[i][j];
		    v[j] = (*this)[j][j] * d;
		    sum -= v[j] * d;
		}

		if ( sum == 0.0f ) {
			return false;
		}

		(*this)[i][i] = sum;
		d = 1.0f / sum;

		for ( j = i + 1; j < numRows; j++ ) {
		    sum = (*this)[j][i];
			for ( k = 0; k < i; k++ ) {
				sum -= (*this)[j][k] * v[k];
			}
		    (*this)[j][i] = sum * d;
		}
	}

	return true;
}

/*
============
CMatXD::ldlt_UpdateRankOne

  Updates the in-place LDL' factorization to obtain the factors for the matrix: LDL' + alpha * v * v'
  If offset > 0 only the lower right corner starting at (offset, offset) is updated.
============
*/
bool CMatXD::ldlt_UpdateRankOne( const CVecXD &v, float alpha, int offset ) {
	int i, j;
	float *y;
	double diag, newDiag, beta, p, d;

	SMF_ASSERT( numRows == numColumns );
	SMF_ASSERT( v.getSize() >= numRows );
	SMF_ASSERT( offset >= 0 && offset < numRows );

	y = (float *) _allocafloat16( v.getSize() * sizeof( float ) );
	memcpy( y, v.toFloatPtr(), v.getSize() * sizeof( float ) );

	for ( i = offset; i < numColumns; i++ ) {
		p = y[i];
		diag = (*this)[i][i];
		(*this)[i][i] = newDiag = diag + alpha * p * p;

		if ( newDiag == 0.0f ) {
			return false;
		}

		alpha /= newDiag;
		beta = p * alpha;
		alpha *= diag;

		for ( j = i+1; j < numRows; j++ ) {

			d = (*this)[j][i];

			y[j] -= p * d;
			d += beta * y[j];

			(*this)[j][i] = d;
		}
	}

	return true;
}

/*
============
CMatXD::ldlt_UpdateRowColumn

  Updates the in-place LDL' factorization to obtain the factors for the matrix:

         [ 0  a  0 ]
  LDL' + [ a  b  c ]
         [ 0  c  0 ]

  where: a = v[0,r-1], b = v[r], c = v[r+1,numRows-1]
============
*/
bool CMatXD::ldlt_UpdateRowColumn( const CVecXD &v, int r ) {
	int i, j;
	double sum;
	float *original, *y;
	CVecXD addSub;

	SMF_ASSERT( numRows == numColumns );
	SMF_ASSERT( v.getSize() >= numRows );
	SMF_ASSERT( r >= 0 && r < numRows );

	addSub.setData( numColumns, (float *) _allocafloat16( numColumns * sizeof( float ) ) );

	if ( r == 0 ) {

		if ( numColumns == 1 ) {
			(*this)[0][0] += v[0];
			return true;
		}
		for ( i = 0; i < numColumns; i++ ) {
			addSub[i] = v[i];
		}

	} else {

		original = (float *) _allocafloat16( numColumns * sizeof( float ) );
		y = (float *) _allocafloat16( numColumns * sizeof( float ) );

		// calculate original row/column of matrix
		for ( i = 0; i < r; i++ ) {
			y[i] = (*this)[r][i] * (*this)[i][i];
		}
		for ( i = 0; i < numColumns; i++ ) {
			if ( i < r ) {
				sum = (*this)[i][i] * (*this)[r][i];
			} else if ( i == r ) {
				sum = (*this)[r][r];
			} else {
				sum = (*this)[r][r] * (*this)[i][r];
			}
			for ( j = 0; j < i && j < r; j++ ) {
				sum += (*this)[i][j] * y[j];
			}
			original[i] = sum;
		}

		// solve for y in L * y = original + v
		for ( i = 0; i < r; i++ ) {
			sum = original[i] + v[i];
			for ( j = 0; j < i; j++ ) {
				sum -= (*this)[i][j] * y[j];
			}
			y[i] = sum;
		}

		// calculate new row of L
		for ( i = 0; i < r; i++ ) {
			(*this)[r][i] = y[i] / (*this)[i][i];
		}

		// if the last row/column of the matrix is updated
		if ( r == numColumns - 1 ) {
			// only calculate new diagonal
			sum = original[r] + v[r];
			for ( j = 0; j < r; j++ ) {
				sum -= (*this)[r][j] * y[j];
			}
			if ( sum == 0.0f ) {
				return false;
			}
			(*this)[r][r] = sum;
			return true;
		}

		// calculate the row/column to be added to the lower right sub matrix starting at (r, r)
		for ( i = 0; i < r; i++ ) {
			y[i] = (*this)[r][i] * (*this)[i][i];
		}
		for ( i = r; i < numColumns; i++ ) {
			if ( i == r ) {
				sum = (*this)[r][r];
			} else {
				sum = (*this)[r][r] * (*this)[i][r];
			}
			for ( j = 0; j < r; j++ ) {
				sum += (*this)[i][j] * y[j];
			}
			addSub[i] = v[i] - ( sum - original[i] );
		}
	}

	// add row/column to the lower right sub matrix starting at (r, r)

#if 0

	CVecXD v1, v2;
	double d;

	v1.setData( numColumns, (float *) _alloca16( numColumns * sizeof( float ) ) );
	v2.setData( numColumns, (float *) _alloca16( numColumns * sizeof( float ) ) );

	d = CMath::SQRT_1OVER2;
	v1[r] = ( 0.5f * addSub[r] + 1.0f ) * d;
	v2[r] = ( 0.5f * addSub[r] - 1.0f ) * d;
	for ( i = r+1; i < numColumns; i++ ) {
		v1[i] = v2[i] = addSub[i] * d;
	}

	// update
	if ( !ldlt_UpdateRankOne( v1, 1.0f, r ) ) {
		return false;
	}
	// downdate
	if ( !ldlt_UpdateRankOne( v2, -1.0f, r ) ) {
		return false;
	}

#else

	float *v1, *v2;
	double d, diag, newDiag, p1, p2, alpha1, alpha2, beta1, beta2;

	v1 = (float *) _allocafloat16( numColumns * sizeof( float ) );
	v2 = (float *) _allocafloat16( numColumns * sizeof( float ) );

	d = CMath::SQRT_1OVER2;
	v1[r] = ( 0.5f * addSub[r] + 1.0f ) * d;
	v2[r] = ( 0.5f * addSub[r] - 1.0f ) * d;
	for ( i = r+1; i < numColumns; i++ ) {
		v1[i] = v2[i] = addSub[i] * d;
	}

	alpha1 = 1.0f;
	alpha2 = -1.0f;

	// simultaneous update/downdate of the sub matrix starting at (r, r)
	for ( i = r; i < numColumns; i++ ) {

		diag = (*this)[i][i];
		p1 = v1[i];
		newDiag = diag + alpha1 * p1 * p1;

		if ( newDiag == 0.0f ) {
			return false;
		}

		alpha1 /= newDiag;
		beta1 = p1 * alpha1;
		alpha1 *= diag;

		diag = newDiag;
		p2 = v2[i];
		newDiag = diag + alpha2 * p2 * p2;

		if ( newDiag == 0.0f ) {
			return false;
		}

		alpha2 /= newDiag;
		beta2 = p2 * alpha2;
		alpha2 *= diag;

		(*this)[i][i] = newDiag;

		for ( j = i+1; j < numRows; j++ ) {

			d = (*this)[j][i];

			v1[j] -= p1 * d;
			d += beta1 * v1[j];

			v2[j] -= p2 * d;
			d += beta2 * v2[j];

			(*this)[j][i] = d;
		}
	}

#endif

	return true;
}

/*
============
CMatXD::ldlt_UpdateIncrement

  Updates the in-place LDL' factorization to obtain the factors for the matrix:

  [ A  a ]
  [ a  b ]

  where: a = v[0,numRows-1], b = v[numRows]
============
*/
bool CMatXD::ldlt_UpdateIncrement( const CVecXD &v ) {
	int i, j;
	float *x;
	double sum, d;

	SMF_ASSERT( numRows == numColumns );
	SMF_ASSERT( v.getSize() >= numRows+1 );

	changeSize( numRows+1, numColumns+1, false );

	x = (float *) _allocafloat16( numRows * sizeof( float ) );

	// solve for x in L * x = v
	for ( i = 0; i < numRows - 1; i++ ) {
		sum = v[i];
		for ( j = 0; j < i; j++ ) {
			sum -= (*this)[i][j] * x[j];
		}
		x[i] = sum;
	}

	// calculate new row of L and calculate the diagonal entry
	sum = v[numRows - 1];
	for ( i = 0; i < numRows - 1; i++ ) {
		(*this)[numRows - 1][i] = d = x[i] / (*this)[i][i];
		sum -= d * x[i];
	}

	if ( sum == 0.0f ) {
		return false;
	}

	// store the diagonal entry
	(*this)[numRows - 1][numRows - 1] = sum;

	return true;
}

/*
============
CMatXD::ldlt_UpdateDecrement

  Updates the in-place LDL' factorization to obtain the factors for the matrix with row r and column r removed.
  v should store the row of the original matrix.
============
*/
bool CMatXD::ldlt_UpdateDecrement( const CVecXD &v, int r ) {
	CVecXD v1;

	SMF_ASSERT( numRows == numColumns );
	SMF_ASSERT( v.getSize() >= numRows );
	SMF_ASSERT( r >= 0 && r < numRows );

	v1.setData( numRows, VECX_ALLOCAFLOAT( numRows ) );

	// update the row and column to identity
	v1 = -v;
	v1[r] += 1.0f;

	// NOTE:	msvc compiler bug: the this pointer stored in edi is expected to stay
	//			untouched when calling ldlt_UpdateRowColumn in the if statement
#if 0
	if ( !ldlt_UpdateRowColumn( v1, r ) ) {
#else
	bool ret = ldlt_UpdateRowColumn( v1, r );
	if ( !ret ) {
#endif
		return false;
	}

	// physically remove the row and column
	update_Decrement( r );

	return true;
}

/*
============
CMatXD::ldlt_solve

  solve Ax = b with A factored in-place as: LDL'
============
*/
void CMatXD::ldlt_solve( CVecXD &x, const CVecXD &b ) const {
	int i, j;
	double sum;

	SMF_ASSERT( numRows == numColumns );
	SMF_ASSERT( x.getSize() >= numRows && b.getSize() >= numRows );

	// solve L
	for ( i = 0; i < numRows; i++ ) {
		sum = b[i];
		for ( j = 0; j < i; j++ ) {
			sum -= (*this)[i][j] * x[j];
		}
		x[i] = sum;
	}

	// solve D
	for ( i = 0; i < numRows; i++ ) {
		x[i] /= (*this)[i][i];
	}

	// solve Lt
	for ( i = numRows - 2; i >= 0; i-- ) {
		sum = x[i];
		for ( j = i + 1; j < numRows; j++ ) {
			sum -= (*this)[j][i] * x[j];
		}
		x[i] = sum;
	}
}

/*
============
CMatXD::ldlt_Inverse

  Calculates the inverse of the matrix which is factored in-place as: LDL'
============
*/
void CMatXD::ldlt_Inverse( CMatXD &inv ) const {
	int i, j;
	CVecXD x, b;

	SMF_ASSERT( numRows == numColumns );

	x.setData( numRows, VECX_ALLOCAFLOAT( numRows ) );
	b.setData( numRows, VECX_ALLOCAFLOAT( numRows ) );
	b.toZero();
	inv.setSize( numRows, numColumns );

	for ( i = 0; i < numRows; i++ ) {

		b[i] = 1.0f;
		ldlt_solve( x, b );
		for ( j = 0; j < numRows; j++ ) {
			inv[j][i] = x[j];
		}
		b[i] = 0.0f;
	}
}

/*
============
CMatXD::ldlt_UnpackFactors

  Unpacks the in-place LDL' factorization.
============
*/
void CMatXD::ldlt_UnpackFactors( CMatXD &L, CMatXD &D ) const {
	int i, j;

	L.zero( numRows, numColumns );
	D.zero( numRows, numColumns );
	for ( i = 0; i < numRows; i++ ) {
		for ( j = 0; j < i; j++ ) {
			L[i][j] = (*this)[i][j];
		}
		L[i][i] = 1.0f;
		D[i][i] = (*this)[i][i];
	}
}

/*
============
CMatXD::ldlt_MultiplyFactors

  Multiplies the factors of the in-place LDL' factorization to form the original matrix.
============
*/
void CMatXD::ldlt_MultiplyFactors( CMatXD &m ) const {
	int r, i, j;
	float *v;
	double sum;

	v = (float *) _allocafloat16( numRows * sizeof( float ) );
	m.setSize( numRows, numColumns );

	for ( r = 0; r < numRows; r++ ) {

		// calculate row of matrix
		for ( i = 0; i < r; i++ ) {
			v[i] = (*this)[r][i] * (*this)[i][i];
		}
		for ( i = 0; i < numColumns; i++ ) {
			if ( i < r ) {
				sum = (*this)[i][i] * (*this)[r][i];
			} else if ( i == r ) {
				sum = (*this)[r][r];
			} else {
				sum = (*this)[r][r] * (*this)[i][r];
			}
			for ( j = 0; j < i && j < r; j++ ) {
				sum += (*this)[i][j] * v[j];
			}
			m[r][i] = sum;
		}
	}
}

/*
============
CMatXD::triDiagonal_ClearTriangles
============
*/
void CMatXD::triDiagonal_ClearTriangles() {
	int i, j;

	SMF_ASSERT( numRows == numColumns );
	for ( i = 0; i < numRows-2; i++ ) {
		for ( j = i+2; j < numColumns; j++ ) {
			(*this)[i][j] = 0.0f;
			(*this)[j][i] = 0.0f;
		}
	}
}

/*
============
CMatXD::triDiagonal_solve

  solve Ax = b with A being tridiagonal.
============
*/
bool CMatXD::triDiagonal_solve( CVecXD &x, const CVecXD &b ) const {
	int i;
	float d;
	CVecXD tmp;

	SMF_ASSERT( numRows == numColumns );
	SMF_ASSERT( x.getSize() >= numRows && b.getSize() >= numRows );

	tmp.setData( numRows, VECX_ALLOCAFLOAT( numRows ) );

	d = (*this)[0][0];
	if ( d == 0.0f ) {
		return false;
	}
	d = 1.0f / d;
	x[0] = b[0] * d;
	for ( i = 1; i < numRows; i++ ) {
		tmp[i] = (*this)[i-1][i] * d;
		d = (*this)[i][i] - (*this)[i][i-1] * tmp[i];
		if ( d == 0.0f ) {
			return false;
		}
		d = 1.0f / d;
		x[i] = ( b[i] - (*this)[i][i-1] * x[i-1] ) * d;
	}
	for ( i = numRows - 2; i >= 0; i-- ) {
		x[i] -= tmp[i+1] * x[i+1];
	}
	return true;
}

/*
============
CMatXD::triDiagonal_Inverse

  Calculates the inverse of a tri-diagonal matrix.
============
*/
void CMatXD::triDiagonal_Inverse( CMatXD &inv ) const {
	int i, j;
	CVecXD x, b;

	SMF_ASSERT( numRows == numColumns );

	x.setData( numRows, VECX_ALLOCAFLOAT( numRows ) );
	b.setData( numRows, VECX_ALLOCAFLOAT( numRows ) );
	b.toZero();
	inv.setSize( numRows, numColumns );

	for ( i = 0; i < numRows; i++ ) {

		b[i] = 1.0f;
		triDiagonal_solve( x, b );
		for ( j = 0; j < numRows; j++ ) {
			inv[j][i] = x[j];
		}
		b[i] = 0.0f;
	}
}

/*
============
CMatXD::HouseholderReduction

  Householder reduction to symmetric tri-diagonal form.
  The original matrix is replaced by an orthogonal matrix effecting the accumulated householder transformations.
  The diagonal elements of the diagonal matrix are stored in diag.
  The off-diagonal elements of the diagonal matrix are stored in subd.
  The initial matrix has to be symmetric.
============
*/
void CMatXD::HouseholderReduction( CVecXD &diag, CVecXD &subd ) {
	int i0, i1, i2, i3;
	float h, f, g, invH, halfFdivH, scale, invScale, sum;

	SMF_ASSERT( numRows == numColumns );

	diag.setSize( numRows );
	subd.setSize( numRows );

	for ( i0 = numRows-1, i3 = numRows-2; i0 >= 1; i0--, i3-- ) {
		h = 0.0f;
		scale = 0.0f;

		if ( i3 > 0 ) {
			for ( i2 = 0; i2 <= i3; i2++ ) {
				scale += CMath::fabs( (*this)[i0][i2] );
			}
			if ( scale == 0 ) {
				subd[i0] = (*this)[i0][i3];
			} else {
				invScale = 1.0f / scale;
				for (i2 = 0; i2 <= i3; i2++)
				{
					(*this)[i0][i2] *= invScale;
					h += (*this)[i0][i2] * (*this)[i0][i2];
				}
				f = (*this)[i0][i3];
				g = CMath::sqrt( h );
				if ( f > 0.0f ) {
					g = -g;
				}
				subd[i0] = scale * g;
				h -= f * g;
				(*this)[i0][i3] = f - g;
				f = 0.0f;
				invH = 1.0f / h;
				for (i1 = 0; i1 <= i3; i1++) {
					(*this)[i1][i0] = (*this)[i0][i1] * invH;
					g = 0.0f;
					for (i2 = 0; i2 <= i1; i2++) {
						g += (*this)[i1][i2] * (*this)[i0][i2];
					}
					for (i2 = i1+1; i2 <= i3; i2++) {
						g += (*this)[i2][i1] * (*this)[i0][i2];
					}
					subd[i1] = g * invH;
					f += subd[i1] * (*this)[i0][i1];
				}
				halfFdivH = 0.5f * f * invH;
				for ( i1 = 0; i1 <= i3; i1++ ) {
					f = (*this)[i0][i1];
					g = subd[i1] - halfFdivH * f;
					subd[i1] = g;
					for ( i2 = 0; i2 <= i1; i2++ ) {
						(*this)[i1][i2] -= f * subd[i2] + g * (*this)[i0][i2];
					}
				}
            }
		} else {
			subd[i0] = (*this)[i0][i3];
		}

		diag[i0] = h;
	}

	diag[0] = 0.0f;
	subd[0] = 0.0f;
	for ( i0 = 0, i3 = -1; i0 <= numRows-1; i0++, i3++ ) {
		if ( diag[i0] ) {
			for ( i1 = 0; i1 <= i3; i1++ ) {
				sum = 0.0f;
				for (i2 = 0; i2 <= i3; i2++) {
					sum += (*this)[i0][i2] * (*this)[i2][i1];
				}
				for ( i2 = 0; i2 <= i3; i2++ ) {
					(*this)[i2][i1] -= sum * (*this)[i2][i0];
				}
			}
		}
		diag[i0] = (*this)[i0][i0];
		(*this)[i0][i0] = 1.0f;
		for ( i1 = 0; i1 <= i3; i1++ ) {
			(*this)[i1][i0] = 0.0f;
			(*this)[i0][i1] = 0.0f;
		}
	}

	// re-order
	for ( i0 = 1, i3 = 0; i0 < numRows; i0++, i3++ ) {
		subd[i3] = subd[i0];
	}
	subd[numRows-1] = 0.0f;
}

/*
============
CMatXD::QL

  QL algorithm with implicit shifts to determine the eigenvalues and eigenvectors of a symmetric tri-diagonal matrix.
  diag contains the diagonal elements of the symmetric tri-diagonal matrix on input and is overwritten with the eigenvalues.
  subd contains the off-diagonal elements of the symmetric tri-diagonal matrix and is destroyed.
  This matrix has to be either the identity matrix to determine the eigenvectors for a symmetric tri-diagonal matrix,
  or the matrix returned by the Householder reduction to determine the eigenvalues for the original symmetric matrix.
============
*/
bool CMatXD::QL( CVecXD &diag, CVecXD &subd ) {
    const int maxIter = 32;
	int i0, i1, i2, i3;
	float a, b, f, g, r, p, s, c;

	SMF_ASSERT( numRows == numColumns );

	for ( i0 = 0; i0 < numRows; i0++ ) {
		for ( i1 = 0; i1 < maxIter; i1++ ) {
			for ( i2 = i0; i2 <= numRows - 2; i2++ ) {
				a = CMath::fabs( diag[i2] ) + CMath::fabs( diag[i2+1] );
				if ( CMath::fabs( subd[i2] ) + a == a ) {
					break;
				}
			}
			if ( i2 == i0 ) {
				break;
			}

			g = ( diag[i0+1] - diag[i0] ) / ( 2.0f * subd[i0] );
			r = CMath::sqrt( g * g + 1.0f );
			if ( g < 0.0f ) {
				g = diag[i2] - diag[i0] + subd[i0] / ( g - r );
			} else {
				g = diag[i2] - diag[i0] + subd[i0] / ( g + r );
			}
			s = 1.0f;
			c = 1.0f;
			p = 0.0f;
			for ( i3 = i2 - 1; i3 >= i0; i3-- ) {
				f = s * subd[i3];
				b = c * subd[i3];
				if ( CMath::fabs( f ) >= CMath::fabs( g ) ) {
					c = g / f;
					r = CMath::sqrt( c * c + 1.0f );
					subd[i3+1] = f * r;
					s = 1.0f / r;
					c *= s;
				} else {
					s = f / g;
					r = CMath::sqrt( s * s + 1.0f );
					subd[i3+1] = g * r;
					c = 1.0f / r;
					s *= c;
				}
				g = diag[i3+1] - p;
				r = ( diag[i3] - g ) * s + 2.0f * b * c;
				p = s * r;
				diag[i3+1] = g + p;
				g = c * r - b;

				for ( int i4 = 0; i4 < numRows; i4++ ) {
					f = (*this)[i4][i3+1];
					(*this)[i4][i3+1] = s * (*this)[i4][i3] + c * f;
					(*this)[i4][i3] = c * (*this)[i4][i3] - s * f;
				}
			}
			diag[i0] -= p;
			subd[i0] = g;
			subd[i2] = 0.0f;
		}
		if ( i1 == maxIter ) {
			return false;
		}
	}
	return true;
}

/*
============
CMatXD::eigen_solveSymmetricTriDiagonal

  Determine eigen values and eigen vectors for a symmetric tri-diagonal matrix.
  The eigen values are stored in 'eigenValues'.
  Column i of the original matrix will store the eigen vector corresponding to the eigenValues[i].
  The initial matrix has to be symmetric tri-diagonal.
============
*/
bool CMatXD::eigen_solveSymmetricTriDiagonal( CVecXD &eigenValues ) {
	int i;
	CVecXD subd;

	SMF_ASSERT( numRows == numColumns );

	subd.setData( numRows, VECX_ALLOCAFLOAT( numRows ) );
	eigenValues.setSize( numRows );

	for ( i = 0; i < numRows-1; i++ ) {
		eigenValues[i] = (*this)[i][i];
		subd[i] = (*this)[i+1][i];
	}
	eigenValues[numRows-1] = (*this)[numRows-1][numRows-1];

	toIdentity();

	return QL( eigenValues, subd );
}

/*
============
CMatXD::eigen_solveSymmetric

  Determine eigen values and eigen vectors for a symmetric matrix.
  The eigen values are stored in 'eigenValues'.
  Column i of the original matrix will store the eigen vector corresponding to the eigenValues[i].
  The initial matrix has to be symmetric.
============
*/
bool CMatXD::eigen_solveSymmetric( CVecXD &eigenValues ) {
	CVecXD subd;

	SMF_ASSERT( numRows == numColumns );

	subd.setData( numRows, VECX_ALLOCAFLOAT( numRows ) );
	eigenValues.setSize( numRows );

	HouseholderReduction( eigenValues, subd );
	return QL( eigenValues, subd );
}

/*
============
CMatXD::HessenbergReduction

  Reduction to Hessenberg form.
============
*/
void CMatXD::HessenbergReduction( CMatXD &H ) {
	int i, j, m;
	int low = 0;
	int high = numRows - 1;
	float scale, f, g, h;
	CVecXD v;

	v.setData( numRows, VECX_ALLOCAFLOAT( numRows ) );

	for ( m = low + 1; m <= high - 1; m++ ) {

		scale = 0.0f;
		for ( i = m; i <= high; i++ ) {
			scale = scale + CMath::fabs( H[i][m-1] );
		}
		if ( scale != 0.0f ) {

			// compute Householder transformation.
			h = 0.0f;
			for ( i = high; i >= m; i-- ) {
				v[i] = H[i][m-1] / scale;
				h += v[i] * v[i];
			}
			g = CMath::sqrt( h );
			if ( v[m] > 0.0f ) {
				g = -g;
			}
			h = h - v[m] * g;
			v[m] = v[m] - g;

			// apply Householder similarity transformation
			// H = (I-u*u'/h)*H*(I-u*u')/h)
			for ( j = m; j < numRows; j++) {
				f = 0.0f;
				for ( i = high; i >= m; i-- ) {
					f += v[i] * H[i][j];
				}
				f = f / h;
				for ( i = m; i <= high; i++ ) {
					H[i][j] -= f * v[i];
				}
			}

			for ( i = 0; i <= high; i++ ) {
				f = 0.0f;
				for ( j = high; j >= m; j-- ) {
					f += v[j] * H[i][j];
				}
				f = f / h;
				for ( j = m; j <= high; j++ ) {
					H[i][j] -= f * v[j];
				}
			}
			v[m] = scale * v[m];
			H[m][m-1] = scale * g;
		}
	}

	// accumulate transformations
	toIdentity();
	for ( int m = high - 1; m >= low + 1; m-- ) {
		if ( H[m][m-1] != 0.0f ) {
			for ( i = m + 1; i <= high; i++ ) {
				v[i] = H[i][m-1];
			}
			for ( j = m; j <= high; j++ ) {
				g = 0.0f;
				for ( i = m; i <= high; i++ ) {
					g += v[i] * (*this)[i][j];
				}
				// float division to avoid possible underflow
				g = ( g / v[m] ) / H[m][m-1];
				for ( i = m; i <= high; i++ ) {
					(*this)[i][j] += g * v[i];
				}
			}
		}
	}
}

/*
============
CMatXD::ComplexDivision

  Complex scalar division.
============
*/
void CMatXD::ComplexDivision( float xr, float xi, float yr, float yi, float &cdivr, float &cdivi ) {
	float r, d;
	if ( CMath::fabs( yr ) > CMath::fabs( yi ) ) {
		r = yi / yr;
		d = yr + r * yi;
		cdivr = ( xr + r * xi ) / d;
		cdivi = ( xi - r * xr ) / d;
	} else {
		r = yr / yi;
		d = yi + r * yr;
		cdivr = ( r * xr + xi ) / d;
		cdivi = ( r * xi - xr ) / d;
	}
}

/*
============
CMatXD::HessenbergToRealSchur

  Reduction from Hessenberg to real Schur form.
============
*/
bool CMatXD::HessenbergToRealSchur( CMatXD &H, CVecXD &realEigenValues, CVecXD &imaginaryEigenValues ) {
	int i, j, k;
	int n = numRows - 1;
	int low = 0;
	int high = numRows - 1;
	float eps = 2e-16f, exshift = 0.0f;
	float p = 0.0f, q = 0.0f, r = 0.0f, s = 0.0f, z = 0.0f, t, w, x, y;

	// store roots isolated by balanc and compute matrix norm
	float norm = 0.0f;
	for ( i = 0; i < numRows; i++ ) {
		if ( i < low || i > high ) {
			realEigenValues[i] = H[i][i];
			imaginaryEigenValues[i] = 0.0f;
		}
		for ( j = MAX( i - 1, 0 ); j < numRows; j++ ) {
			norm = norm + CMath::fabs( H[i][j] );
		}
	}

	int iter = 0;
	while( n >= low ) {

		// look for single small sub-diagonal element
		int l = n;
		while ( l > low ) {
			s = CMath::fabs( H[l-1][l-1] ) + CMath::fabs( H[l][l] );
			if ( s == 0.0f ) {
				s = norm;
			}
			if ( CMath::fabs( H[l][l-1] ) < eps * s ) {
				break;
			}
			l--;
		}

		// check for convergence
		if ( l == n ) {			// one root found
			H[n][n] = H[n][n] + exshift;
			realEigenValues[n] = H[n][n];
			imaginaryEigenValues[n] = 0.0f;
			n--;
			iter = 0;
		} else if ( l == n-1 ) {	// two roots found
			w = H[n][n-1] * H[n-1][n];
			p = ( H[n-1][n-1] - H[n][n] ) / 2.0f;
			q = p * p + w;
			z = CMath::sqrt( CMath::fabs( q ) );
			H[n][n] = H[n][n] + exshift;
			H[n-1][n-1] = H[n-1][n-1] + exshift;
			x = H[n][n];

			if ( q >= 0.0f ) {		// real pair
				if ( p >= 0.0f ) {
					z = p + z;
				} else {
					z = p - z;
				}
				realEigenValues[n-1] = x + z;
				realEigenValues[n] = realEigenValues[n-1];
				if ( z != 0.0f ) {
					realEigenValues[n] = x - w / z;
				}
				imaginaryEigenValues[n-1] = 0.0f;
				imaginaryEigenValues[n] = 0.0f;
				x = H[n][n-1];
				s = CMath::fabs( x ) + CMath::fabs( z );
				p = x / s;
				q = z / s;
				r = CMath::sqrt( p * p + q * q );
				p = p / r;
				q = q / r;

				// modify row
				for ( j = n-1; j < numRows; j++ ) {
					z = H[n-1][j];
					H[n-1][j] = q * z + p * H[n][j];
					H[n][j] = q * H[n][j] - p * z;
				}

				// modify column
				for ( i = 0; i <= n; i++ ) {
					z = H[i][n-1];
					H[i][n-1] = q * z + p * H[i][n];
					H[i][n] = q * H[i][n] - p * z;
				}

				// accumulate transformations
				for ( i = low; i <= high; i++ ) {
					z = (*this)[i][n-1];
					(*this)[i][n-1] = q * z + p * (*this)[i][n];
					(*this)[i][n] = q * (*this)[i][n] - p * z;
				}
			} else {		// complex pair
				realEigenValues[n-1] = x + p;
				realEigenValues[n] = x + p;
				imaginaryEigenValues[n-1] = z;
				imaginaryEigenValues[n] = -z;
			}
			n = n - 2;
			iter = 0;

		} else {	// no convergence yet

			// form shift
			x = H[n][n];
			y = 0.0f;
			w = 0.0f;
			if ( l < n ) {
				y = H[n-1][n-1];
				w = H[n][n-1] * H[n-1][n];
			}

			// Wilkinson's original ad hoc shift
			if ( iter == 10 ) {
				exshift += x;
				for ( i = low; i <= n; i++ ) {
					H[i][i] -= x;
				}
				s = CMath::fabs( H[n][n-1] ) + CMath::fabs( H[n-1][n-2] );
				x = y = 0.75f * s;
				w = -0.4375f * s * s;
			}

			// new ad hoc shift
			if ( iter == 30 ) {
				s = ( y - x ) / 2.0f;
				s = s * s + w;
				if ( s > 0 ) {
					s = CMath::sqrt( s );
					if ( y < x ) {
						s = -s;
					}
					s = x - w / ( ( y - x ) / 2.0f + s );
					for ( i = low; i <= n; i++ ) {
						H[i][i] -= s;
					}
					exshift += s;
					x = y = w = 0.964f;
				}
			}

			iter = iter + 1;

			// look for two consecutive small sub-diagonal elements
			int m;
			for( m = n-2; m >= l; m-- ) {
				z = H[m][m];
				r = x - z;
				s = y - z;
				p = ( r * s - w ) / H[m+1][m] + H[m][m+1];
				q = H[m+1][m+1] - z - r - s;
				r = H[m+2][m+1];
				s = CMath::fabs( p ) + CMath::fabs( q ) + CMath::fabs( r );
				p = p / s;
				q = q / s;
				r = r / s;
				if ( m == l ) {
					break;
				}
				if ( CMath::fabs( H[m][m-1] ) * ( CMath::fabs( q ) + CMath::fabs( r ) ) <
						eps * ( CMath::fabs( p ) * ( CMath::fabs( H[m-1][m-1] ) + CMath::fabs( z ) + CMath::fabs( H[m+1][m+1] ) ) ) ) {
					break;
				}
			}

			for ( i = m+2; i <= n; i++ ) {
				H[i][i-2] = 0.0f;
				if ( i > m+2 ) {
					H[i][i-3] = 0.0f;
				}
			}

			// double QR step involving rows l:n and columns m:n
			for ( k = m; k <= n-1; k++ ) {
				bool notlast = ( k != n-1 );
				if ( k != m ) {
					p = H[k][k-1];
					q = H[k+1][k-1];
					r = ( notlast ? H[k+2][k-1] : 0.0f );
					x = CMath::fabs( p ) + CMath::fabs( q ) + CMath::fabs( r );
					if ( x != 0.0f ) {
						p = p / x;
						q = q / x;
						r = r / x;
					}
				}
				if ( x == 0.0f ) {
					break;
				}
				s = CMath::sqrt( p * p + q * q + r * r );
				if ( p < 0.0f ) {
					s = -s;
				}
				if ( s != 0.0f ) {
					if ( k != m ) {
						H[k][k-1] = -s * x;
					} else if ( l != m ) {
						H[k][k-1] = -H[k][k-1];
					}
					p = p + s;
					x = p / s;
					y = q / s;
					z = r / s;
					q = q / p;
					r = r / p;

					// modify row
					for ( j = k; j < numRows; j++ ) {
						p = H[k][j] + q * H[k+1][j];
						if ( notlast ) {
							p = p + r * H[k+2][j];
							H[k+2][j] = H[k+2][j] - p * z;
						}
						H[k][j] = H[k][j] - p * x;
						H[k+1][j] = H[k+1][j] - p * y;
					}

					// modify column
					for ( i = 0; i <= MIN( n, k + 3 ); i++ ) {
						p = x * H[i][k] + y * H[i][k+1];
						if ( notlast ) {
							p = p + z * H[i][k+2];
							H[i][k+2] = H[i][k+2] - p * r;
						}
						H[i][k] = H[i][k] - p;
						H[i][k+1] = H[i][k+1] - p * q;
					}

					// accumulate transformations
					for ( i = low; i <= high; i++ ) {
						p = x * (*this)[i][k] + y * (*this)[i][k+1];
						if ( notlast ) {
							p = p + z * (*this)[i][k+2];
							(*this)[i][k+2] = (*this)[i][k+2] - p * r;
						}
						(*this)[i][k] = (*this)[i][k] - p;
						(*this)[i][k+1] = (*this)[i][k+1] - p * q;
					}
				}
			}
		}
	}

	// backsubstitute to find vectors of upper triangular form
	if ( norm == 0.0f ) {
		return false;
	}

	for ( n = numRows-1; n >= 0; n-- ) {
		p = realEigenValues[n];
		q = imaginaryEigenValues[n];

		if ( q == 0.0f ) {		// real vector
			int l = n;
			H[n][n] = 1.0f;
			for ( i = n-1; i >= 0; i-- ) {
				w = H[i][i] - p;
				r = 0.0f;
				for ( j = l; j <= n; j++ ) {
					r = r + H[i][j] * H[j][n];
				}
				if ( imaginaryEigenValues[i] < 0.0f ) {
					z = w;
					s = r;
				} else {
					l = i;
					if ( imaginaryEigenValues[i] == 0.0f ) {
						if ( w != 0.0f ) {
							H[i][n] = -r / w;
						} else {
							H[i][n] = -r / ( eps * norm );
						}
					} else {		// solve real equations
						x = H[i][i+1];
						y = H[i+1][i];
						q = ( realEigenValues[i] - p ) * ( realEigenValues[i] - p ) + imaginaryEigenValues[i] * imaginaryEigenValues[i];
						t = ( x * s - z * r ) / q;
						H[i][n] = t;
						if ( CMath::fabs(x) > CMath::fabs( z ) ) {
							H[i+1][n] = ( -r - w * t ) / x;
						} else {
							H[i+1][n] = ( -s - y * t ) / z;
						}
					}

					// overflow control
					t = CMath::fabs(H[i][n]);
					if ( ( eps * t ) * t > 1 ) {
						for ( j = i; j <= n; j++ ) {
							H[j][n] = H[j][n] / t;
						}
					}
				}
			}
		} else if ( q < 0.0f ) {	// complex vector
			int l = n-1;

			// last vector component imaginary so matrix is triangular
			if ( CMath::fabs( H[n][n-1] ) > CMath::fabs( H[n-1][n] ) ) {
				H[n-1][n-1] = q / H[n][n-1];
				H[n-1][n] = -( H[n][n] - p ) / H[n][n-1];
			} else {
				ComplexDivision( 0.0f, -H[n-1][n], H[n-1][n-1]-p, q, H[n-1][n-1], H[n-1][n] );
			}
			H[n][n-1] = 0.0f;
			H[n][n] = 1.0f;
			for ( i = n-2; i >= 0; i-- ) {
				float ra, sa, vr, vi;
				ra = 0.0f;
				sa = 0.0f;
				for ( j = l; j <= n; j++ ) {
					ra = ra + H[i][j] * H[j][n-1];
					sa = sa + H[i][j] * H[j][n];
				}
				w = H[i][i] - p;

				if ( imaginaryEigenValues[i] < 0.0f ) {
					z = w;
					r = ra;
					s = sa;
				} else {
					l = i;
					if ( imaginaryEigenValues[i] == 0.0f ) {
						ComplexDivision( -ra, -sa, w, q, H[i][n-1], H[i][n] );
					} else {
						// solve complex equations
						x = H[i][i+1];
						y = H[i+1][i];
						vr = ( realEigenValues[i] - p ) * ( realEigenValues[i] - p ) + imaginaryEigenValues[i] * imaginaryEigenValues[i] - q * q;
						vi = ( realEigenValues[i] - p ) * 2.0f * q;
						if ( vr == 0.0f && vi == 0.0f ) {
							vr = eps * norm * ( CMath::fabs( w ) + CMath::fabs( q ) + CMath::fabs( x ) + CMath::fabs( y ) + CMath::fabs( z ) );
						}
						ComplexDivision( x * r - z * ra + q * sa, x * s - z * sa - q * ra, vr, vi, H[i][n-1], H[i][n] );
						if ( CMath::fabs( x ) > ( CMath::fabs( z ) + CMath::fabs( q ) ) ) {
							H[i+1][n-1] = ( -ra - w * H[i][n-1] + q * H[i][n] ) / x;
							H[i+1][n] = ( -sa - w * H[i][n] - q * H[i][n-1] ) / x;
						} else {
							ComplexDivision( -r - y * H[i][n-1], -s - y * H[i][n], z, q, H[i+1][n-1], H[i+1][n] );
						}
					}

					// overflow control
					t = MAX( CMath::fabs( H[i][n-1] ), CMath::fabs( H[i][n] ) );
					if ( ( eps * t ) * t > 1 ) {
						for ( j = i; j <= n; j++ ) {
							H[j][n-1] = H[j][n-1] / t;
							H[j][n] = H[j][n] / t;
						}
					}
				}
			}
		}
	}

	// vectors of isolated roots
	for ( i = 0; i < numRows; i++ ) {
		if ( i < low || i > high ) {
			for ( j = i; j < numRows; j++ ) {
				(*this)[i][j] = H[i][j];
			}
		}
	}

	// back transformation to get eigenvectors of original matrix
	for ( j = numRows - 1; j >= low; j-- ) {
		for ( i = low; i <= high; i++ ) {
			z = 0.0f;
			for ( k = low; k <= MIN( j, high ); k++ ) {
				z = z + (*this)[i][k] * H[k][j];
			}
			(*this)[i][j] = z;
		}
	}

	return true;
}

/*
============
CMatXD::eigen_solve

  Determine eigen values and eigen vectors for a square matrix.
  The eigen values are stored in 'realEigenValues' and 'imaginaryEigenValues'.
  Column i of the original matrix will store the eigen vector corresponding to the realEigenValues[i] and imaginaryEigenValues[i].
============
*/
bool CMatXD::eigen_solve( CVecXD &realEigenValues, CVecXD &imaginaryEigenValues ) {
    CMatXD H;

	SMF_ASSERT( numRows == numColumns );

	realEigenValues.setSize( numRows );
	imaginaryEigenValues.setSize( numRows );

	H = *this;

    // reduce to Hessenberg form
    HessenbergReduction( H );

    // reduce Hessenberg to real Schur form
    return HessenbergToRealSchur( H, realEigenValues, imaginaryEigenValues );
}

/*
============
CMatXD::eigen_SortIncreasing
============
*/
void CMatXD::eigen_SortIncreasing( CVecXD &eigenValues ) {
	int i, j, k;
	float min;

	for ( i = 0, j; i <= numRows - 2; i++ ) {
		j = i;
		min = eigenValues[j];
		for ( k = i + 1; k < numRows; k++ ) {
			if ( eigenValues[k] < min ) {
				j = k;
				min = eigenValues[j];
			}
		}
		if ( j != i ) {
			eigenValues.swapElements( i, j );
			swapColumns( i, j );
		}
	}
}

/*
============
CMatXD::eigen_SortDecreasing
============
*/
void CMatXD::eigen_SortDecreasing( CVecXD &eigenValues ) {
	int i, j, k;
	float max;

	for ( i = 0, j; i <= numRows - 2; i++ ) {
		j = i;
		max = eigenValues[j];
		for ( k = i + 1; k < numRows; k++ ) {
			if ( eigenValues[k] > max ) {
				j = k;
				max = eigenValues[j];
			}
		}
		if ( j != i ) {
			eigenValues.swapElements( i, j );
			swapColumns( i, j );
		}
	}
}

/*
============
CMatXD::DeterminantGeneric
============
*/
float CMatXD::DeterminantGeneric() const {
	int *index;
	float det;
	CMatXD tmp;

	index = (int *) _allocafloat16( numRows * sizeof( int ) );
	tmp.setData( numRows, numColumns, MATX_ALLOCAFLOAT( numRows * numColumns ) );
	tmp = *this;

	if ( !tmp.lu_Factor( index, &det ) ) {
		return 0.0f;
	}

	return det;
}

/*
============
CMatXD::InverseSelfGeneric
============
*/
bool CMatXD::InverseSelfGeneric() {
	int i, j, *index;
	CMatXD tmp;
	CVecXD x, b;

	index = (int *) _allocafloat16( numRows * sizeof( int ) );
	tmp.setData( numRows, numColumns, MATX_ALLOCAFLOAT( numRows * numColumns ) );
	tmp = *this;

	if ( !tmp.lu_Factor( index ) ) {
		return false;
	}

	x.setData( numRows, VECX_ALLOCAFLOAT( numRows ) );
	b.setData( numRows, VECX_ALLOCAFLOAT( numRows ) );
	b.toZero();

	for ( i = 0; i < numRows; i++ ) {

		b[i] = 1.0f;
		tmp.lu_solve( x, b, index );
		for ( j = 0; j < numRows; j++ ) {
			(*this)[j][i] = x[j];
		}
		b[i] = 0.0f;
	}
	return true;
}

/*
============
CMatXD::Test
============
*/
void CMatXD::Test() {
	CMatXD original, m1, m2, m3, q1, q2, r1, r2;
	CVecXD v, w, u, c, d;
	int offset, size, *index1, *index2;

	size = 6;
	original.random( size, size, 0 );
	original = original * original.transpose();

	index1 = (int *) _allocafloat16( ( size + 1 ) * sizeof( index1[0] ) );
	index2 = (int *) _allocafloat16( ( size + 1 ) * sizeof( index2[0] ) );

	/*
		CMatXD::lowerTriangularInverse
	*/

	m1 = original;
	m1.clearUpperTriangle();
	m2 = m1;

	m2.inverseSelf();
	m1.lowerTriangularInverse();

	if ( !m1.compare( m2, 1e-4f ) ) {

		Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::lowerTriangularInverse failed" <<endl;
	}

	/*
		CMatXD::upperTriangularInverse
	*/

	m1 = original;
	m1.clearLowerTriangle();
	m2 = m1;

	m2.inverseSelf();
	m1.upperTriangularInverse();

	if ( !m1.compare( m2, 1e-4f ) ) {
		Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::upperTriangularInverse failed" <<endl;
	}

	/*
		CMatXD::inverse_GaussJordan
	*/

	m1 = original;

	m1.inverse_GaussJordan();
	m1 *= original;

	if ( !m1.isIdentity( 1e-4f ) ) {
		Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::inverse_GaussJordan failed" <<endl;
	}

	/*
		CMatXD::inverse_UpdateRankOne
	*/

	m1 = original;
	m2 = original;

	w.random( size, 1 );
	v.random( size, 2 );

	// invert m1
	m1.inverse_GaussJordan();

	// modify and invert m2
	m2.update_RankOne( v, w, 1.0f );
	if ( !m2.inverse_GaussJordan() ) {
		SMF_ASSERT( 0 );
	}

	// update inverse of m1
	m1.inverse_UpdateRankOne( v, w, 1.0f );

	if ( !m1.compare( m2, 1e-4f ) ) {
		Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::inverse_UpdateRankOne failed" <<endl;
	}

	/*
		CMatXD::inverse_UpdateRowColumn
	*/

	for ( offset = 0; offset < size; offset++ ) {
		m1 = original;
		m2 = original;

		v.random( size, 1 );
		w.random( size, 2 );
		w[offset] = 0.0f;

		// invert m1
		m1.inverse_GaussJordan();

		// modify and invert m2
		m2.update_RowColumn( v, w, offset );
		if ( !m2.inverse_GaussJordan() ) {
			SMF_ASSERT( 0 );
		}

		// update inverse of m1
		m1.inverse_UpdateRowColumn( v, w, offset );

		if ( !m1.compare( m2, 1e-3f ) ) {
			Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::inverse_UpdateRowColumn failed" <<endl;
		}
	}

	/*
		CMatXD::inverse_UpdateIncrement
	*/

	m1 = original;
	m2 = original;

	v.random( size + 1, 1 );
	w.random( size + 1, 2 );
	w[size] = 0.0f;

	// invert m1
	m1.inverse_GaussJordan();

	// modify and invert m2
	m2.update_Increment( v, w );
	if ( !m2.inverse_GaussJordan() ) {
		SMF_ASSERT( 0 );
	}

	// update inverse of m1
	m1.inverse_UpdateIncrement( v, w );

	if ( !m1.compare( m2, 1e-4f ) ) {
		Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::inverse_UpdateIncrement failed" <<endl;
	}

	/*
		CMatXD::inverse_UpdateDecrement
	*/

	for ( offset = 0; offset < size; offset++ ) {
		m1 = original;
		m2 = original;

		v.setSize( 6 );
		w.setSize( 6 );
		for ( int i = 0; i < size; i++ ) {
			v[i] = original[i][offset];
			w[i] = original[offset][i];
		}

		// invert m1
		m1.inverse_GaussJordan();

		// modify and invert m2
		m2.update_Decrement( offset );
		if ( !m2.inverse_GaussJordan() ) {
			SMF_ASSERT( 0 );
		}

		// update inverse of m1
		m1.inverse_UpdateDecrement( v, w, offset );

		if ( !m1.compare( m2, 1e-3f ) ) {
			Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::inverse_UpdateDecrement failed" <<endl;
		}
	}

	/*
		CMatXD::lu_Factor
	*/

	m1 = original;

	m1.lu_Factor( NULL );	// no pivoting
	m1.lu_UnpackFactors( m2, m3 );
	m1 = m2 * m3;

	if ( !original.compare( m1, 1e-4f ) ) {
		Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::lu_Factor failed" <<endl;
	}

	/*
		CMatXD::lu_UpdateRankOne
	*/

	m1 = original;
	m2 = original;

	w.random( size, 1 );
	v.random( size, 2 );

	// factor m1
	m1.lu_Factor( index1 );

	// modify and factor m2
	m2.update_RankOne( v, w, 1.0f );
	if ( !m2.lu_Factor( index2 ) ) {
		SMF_ASSERT( 0 );
	}
	m2.lu_MultiplyFactors( m3, index2 );
	m2 = m3;

	// update factored m1
	m1.lu_UpdateRankOne( v, w, 1.0f, index1 );
	m1.lu_MultiplyFactors( m3, index1 );
	m1 = m3;

	if ( !m1.compare( m2, 1e-4f ) ) {
		Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::lu_UpdateRankOne failed" <<endl;
	}

	/*
		CMatXD::lu_UpdateRowColumn
	*/

	for ( offset = 0; offset < size; offset++ ) {
		m1 = original;
		m2 = original;

		v.random( size, 1 );
		w.random( size, 2 );
		w[offset] = 0.0f;

		// factor m1
		m1.lu_Factor( index1 );

		// modify and factor m2
		m2.update_RowColumn( v, w, offset );
		if ( !m2.lu_Factor( index2 ) ) {
			SMF_ASSERT( 0 );
		}
		m2.lu_MultiplyFactors( m3, index2 );
		m2 = m3;

		// update m1
		m1.lu_UpdateRowColumn( v, w, offset, index1  );
		m1.lu_MultiplyFactors( m3, index1 );
		m1 = m3;

		if ( !m1.compare( m2, 1e-3f ) ) {
			Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::lu_UpdateRowColumn failed" <<endl;
		}
	}

	/*
		CMatXD::lu_UpdateIncrement
	*/

	m1 = original;
	m2 = original;

	v.random( size + 1, 1 );
	w.random( size + 1, 2 );
	w[size] = 0.0f;

	// factor m1
	m1.lu_Factor( index1 );

	// modify and factor m2
	m2.update_Increment( v, w );
	if ( !m2.lu_Factor( index2 ) ) {
		SMF_ASSERT( 0 );
	}
	m2.lu_MultiplyFactors( m3, index2 );
	m2 = m3;

	// update factored m1
	m1.lu_UpdateIncrement( v, w, index1 );
	m1.lu_MultiplyFactors( m3, index1 );
	m1 = m3;

	if ( !m1.compare( m2, 1e-4f ) ) {
		Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::lu_UpdateIncrement failed" <<endl;
	}

	/*
		CMatXD::lu_UpdateDecrement
	*/

	for ( offset = 0; offset < size; offset++ ) {
		m1 = original;
		m2 = original;

		v.setSize( 6 );
		w.setSize( 6 );
		for ( int i = 0; i < size; i++ ) {
			v[i] = original[i][offset];
			w[i] = original[offset][i];
		}

		// factor m1
		m1.lu_Factor( index1 );

		// modify and factor m2
		m2.update_Decrement( offset );
		if ( !m2.lu_Factor( index2 ) ) {
			SMF_ASSERT( 0 );
		}
		m2.lu_MultiplyFactors( m3, index2 );
		m2 = m3;

		u.setSize( 6 );
		for ( int i = 0; i < size; i++ ) {
			u[i] = original[index1[offset]][i];
		}

		// update factors of m1
		m1.lu_UpdateDecrement( v, w, u, offset, index1 );
		m1.lu_MultiplyFactors( m3, index1 );
		m1 = m3;

		if ( !m1.compare( m2, 1e-3f ) ) {
			Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::lu_UpdateDecrement failed" <<endl;
		}
	}

	/*
		CMatXD::lu_Inverse
	*/

	m2 = original;

	m2.lu_Factor( NULL );
	m2.lu_Inverse( m1, NULL );
	m1 *= original;

	if ( !m1.isIdentity( 1e-4f ) ) {
		Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::lu_Inverse failed" <<endl;
	}

	/*
		CMatXD::qr_Factor
	*/

	c.setSize( size );
	d.setSize( size );

	m1 = original;

	m1.qr_Factor( c, d );
	m1.qr_UnpackFactors( q1, r1, c, d );
	m1 = q1 * r1;

	if ( !original.compare( m1, 1e-4f ) ) {
		Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::qr_Factor failed" <<endl;
	}

	/*
		CMatXD::qr_UpdateRankOne
	*/

	c.setSize( size );
	d.setSize( size );

	m1 = original;
	m2 = original;

	w.random( size, 0 );
	v = w;

	// factor m1
	m1.qr_Factor( c, d );
	m1.qr_UnpackFactors( q1, r1, c, d );

	// modify and factor m2
	m2.update_RankOne( v, w, 1.0f );
	if ( !m2.qr_Factor( c, d ) ) {
		SMF_ASSERT( 0 );
	}
	m2.qr_UnpackFactors( q2, r2, c, d );
	m2 = q2 * r2;

	// update factored m1
	q1.qr_UpdateRankOne( r1, v, w, 1.0f );
	m1 = q1 * r1;

	if ( !m1.compare( m2, 1e-4f ) ) {
		Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::qr_UpdateRankOne failed" <<endl;
	}

	/*
		CMatXD::qr_UpdateRowColumn
	*/

	for ( offset = 0; offset < size; offset++ ) {
		c.setSize( size );
		d.setSize( size );

		m1 = original;
		m2 = original;

		v.random( size, 1 );
		w.random( size, 2 );
		w[offset] = 0.0f;

		// factor m1
		m1.qr_Factor( c, d );
		m1.qr_UnpackFactors( q1, r1, c, d );

		// modify and factor m2
		m2.update_RowColumn( v, w, offset );
		if ( !m2.qr_Factor( c, d ) ) {
			SMF_ASSERT( 0 );
		}
		m2.qr_UnpackFactors( q2, r2, c, d );
		m2 = q2 * r2;

		// update m1
		q1.qr_UpdateRowColumn( r1, v, w, offset );
		m1 = q1 * r1;

		if ( !m1.compare( m2, 1e-3f ) ) {
			Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::qr_UpdateRowColumn failed" <<endl;
		}
	}

	/*
		CMatXD::qr_UpdateIncrement
	*/

	c.setSize( size+1 );
	d.setSize( size+1 );

	m1 = original;
	m2 = original;

	v.random( size + 1, 1 );
	w.random( size + 1, 2 );
	w[size] = 0.0f;

	// factor m1
	m1.qr_Factor( c, d );
	m1.qr_UnpackFactors( q1, r1, c, d );

	// modify and factor m2
	m2.update_Increment( v, w );
	if ( !m2.qr_Factor( c, d ) ) {
		SMF_ASSERT( 0 );
	}
	m2.qr_UnpackFactors( q2, r2, c, d );
	m2 = q2 * r2;

	// update factored m1
	q1.qr_UpdateIncrement( r1, v, w );
	m1 = q1 * r1;

	if ( !m1.compare( m2, 1e-4f ) ) {
		Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::qr_UpdateIncrement failed" <<endl;
	}

	/*
		CMatXD::qr_UpdateDecrement
	*/

	for ( offset = 0; offset < size; offset++ ) {
		c.setSize( size+1 );
		d.setSize( size+1 );

		m1 = original;
		m2 = original;

		v.setSize( 6 );
		w.setSize( 6 );
		for ( int i = 0; i < size; i++ ) {
			v[i] = original[i][offset];
			w[i] = original[offset][i];
		}

		// factor m1
		m1.qr_Factor( c, d );
		m1.qr_UnpackFactors( q1, r1, c, d );

		// modify and factor m2
		m2.update_Decrement( offset );
		if ( !m2.qr_Factor( c, d ) ) {
			SMF_ASSERT( 0 );
		}
		m2.qr_UnpackFactors( q2, r2, c, d );
		m2 = q2 * r2;

		// update factors of m1
		q1.qr_UpdateDecrement( r1, v, w, offset );
		m1 = q1 * r1;

		if ( !m1.compare( m2, 1e-3f ) ) {
			Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::qr_UpdateDecrement failed" <<endl;
		}
	}

	/*
		CMatXD::qr_Inverse
	*/

	m2 = original;

	m2.qr_Factor( c, d );
	m2.qr_Inverse( m1, c, d );
	m1 *= original;

	if ( !m1.isIdentity( 1e-4f ) ) {
		Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::qr_Inverse failed" <<endl;
	}

	/*
		CMatXD::svd_Factor
	*/

	m1 = original;
	m3.zero( size, size );
	w.zero( size );

	m1.svd_Factor( w, m3 );
	m2.diag( w );
	m3.transposeSelf();
	m1 = m1 * m2 * m3;

	if ( !original.compare( m1, 1e-4f ) ) {
		Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::svd_Factor failed" <<endl;
	}

	/*
		CMatXD::svd_Inverse
	*/

	m2 = original;

	m2.svd_Factor( w, m3 );
	m2.svd_Inverse( m1, w, m3 );
	m1 *= original;

	if ( !m1.isIdentity( 1e-4f ) ) {
		Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::svd_Inverse failed" <<endl;
	}

	/*
		CMatXD::cholesky_Factor
	*/

	m1 = original;

	m1.cholesky_Factor();
	m1.cholesky_MultiplyFactors( m2 );

	if ( !original.compare( m2, 1e-4f ) ) {
		Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::cholesky_Factor failed" <<endl;
	}

	/*
		CMatXD::cholesky_UpdateRankOne
	*/

	m1 = original;
	m2 = original;

	w.random( size, 0 );

	// factor m1
	m1.cholesky_Factor();
	m1.clearUpperTriangle();

	// modify and factor m2
	m2.update_RankOneSymmetric( w, 1.0f );
	if ( !m2.cholesky_Factor() ) {
		SMF_ASSERT( 0 );
	}
	m2.clearUpperTriangle();

	// update factored m1
	m1.cholesky_UpdateRankOne( w, 1.0f, 0 );

	if ( !m1.compare( m2, 1e-4f ) ) {
		Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::cholesky_UpdateRankOne failed" <<endl;
	}

	/*
		CMatXD::cholesky_UpdateRowColumn
	*/

	for ( offset = 0; offset < size; offset++ ) {
		m1 = original;
		m2 = original;

		// factor m1
		m1.cholesky_Factor();
		m1.clearUpperTriangle();

		int pdtable[] = { 1, 0, 1, 0, 0, 0 };
		w.random( size, pdtable[offset] );
		w *= 0.1f;

		// modify and factor m2
		m2.update_RowColumnSymmetric( w, offset );
		if ( !m2.cholesky_Factor() ) {
			SMF_ASSERT( 0 );
		}
		m2.clearUpperTriangle();

		// update m1
		m1.cholesky_UpdateRowColumn( w, offset );

		if ( !m1.compare( m2, 1e-3f ) ) {
			Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::cholesky_UpdateRowColumn failed" <<endl;
		}
	}

	/*
		CMatXD::cholesky_UpdateIncrement
	*/

	m1.random( size + 1, size + 1, 0 );
	m3 = m1 * m1.transpose();

	m1.squareSubMatrix( m3, size );
	m2 = m1;

	w.setSize( size + 1 );
	for ( int i = 0; i < size + 1; i++ ) {
		w[i] = m3[size][i];
	}

	// factor m1
	m1.cholesky_Factor();

	// modify and factor m2
	m2.update_IncrementSymmetric( w );
	if ( !m2.cholesky_Factor() ) {
		SMF_ASSERT( 0 );
	}

	// update factored m1
	m1.cholesky_UpdateIncrement( w );

	m1.clearUpperTriangle();
	m2.clearUpperTriangle();

	if ( !m1.compare( m2, 1e-4f ) ) {
		Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::cholesky_UpdateIncrement failed" <<endl;
	}

	/*
		CMatXD::cholesky_UpdateDecrement
	*/

	for ( offset = 0; offset < size; offset += size - 1 ) {
		m1 = original;
		m2 = original;

		v.setSize( 6 );
		for ( int i = 0; i < size; i++ ) {
			v[i] = original[i][offset];
		}

		// factor m1
		m1.cholesky_Factor();

		// modify and factor m2
		m2.update_Decrement( offset );
		if ( !m2.cholesky_Factor() ) {
			SMF_ASSERT( 0 );
		}

		// update factors of m1
		m1.cholesky_UpdateDecrement( v, offset );

		if ( !m1.compare( m2, 1e-3f ) ) {
			Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::cholesky_UpdateDecrement failed" <<endl;
		}
	}

	/*
		CMatXD::cholesky_Inverse
	*/

	m2 = original;

	m2.cholesky_Factor();
	m2.cholesky_Inverse( m1 );
	m1 *= original;

	if ( !m1.isIdentity( 1e-4f ) ) {
		Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::cholesky_Inverse failed" <<endl;
	}

	/*
		CMatXD::ldlt_Factor
	*/

	m1 = original;

	m1.ldlt_Factor();
	m1.ldlt_MultiplyFactors( m2 );

	if ( !original.compare( m2, 1e-4f ) ) {
		Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::ldlt_Factor failed" <<endl;
	}

	m1.ldlt_UnpackFactors( m2, m3 );
	m2 = m2 * m3 * m2.transpose();

	if ( !original.compare( m2, 1e-4f ) ) {
		Debug::debug(Debug::error,__FUNCTION__)<< "CMatXD::ldlt_Factor failed" <<endl;
	}

	/*
		CMatXD::ldlt_UpdateRankOne
	*/

	m1 = original;
	m2 = original;

	w.random( size, 0 );

	// factor m1
	m1.ldlt_Factor();
	m1.clearUpperTriangle();

	// modify and factor m2
	m2.update_RankOneSymmetric( w, 1.0f );
	if ( !m2.ldlt_Factor() ) {
		SMF_ASSERT( 0 );
	}
	m2.clearUpperTriangle();

	// update factored m1
	m1.ldlt_UpdateRankOne( w, 1.0f, 0 );

	if ( !m1.compare( m2, 1e-4f ) ) {
		Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::ldlt_UpdateRankOne failed" <<endl;
	}

	/*
		CMatXD::ldlt_UpdateRowColumn
	*/

	for ( offset = 0; offset < size; offset++ ) {
		m1 = original;
		m2 = original;

		w.random( size, 0 );

		// factor m1
		m1.ldlt_Factor();
		m1.clearUpperTriangle();

		// modify and factor m2
		m2.update_RowColumnSymmetric( w, offset );
		if ( !m2.ldlt_Factor() ) {
			SMF_ASSERT( 0 );
		}
		m2.clearUpperTriangle();

		// update m1
		m1.ldlt_UpdateRowColumn( w, offset );

		if ( !m1.compare( m2, 1e-3f ) ) {
			Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::ldlt_UpdateRowColumn failed" <<endl;
		}
	}

	/*
		CMatXD::ldlt_UpdateIncrement
	*/

	m1.random( size + 1, size + 1, 0 );
	m3 = m1 * m1.transpose();

	m1.squareSubMatrix( m3, size );
	m2 = m1;

	w.setSize( size + 1 );
	for ( int i = 0; i < size + 1; i++ ) {
		w[i] = m3[size][i];
	}

	// factor m1
	m1.ldlt_Factor();

	// modify and factor m2
	m2.update_IncrementSymmetric( w );
	if ( !m2.ldlt_Factor() ) {
		SMF_ASSERT( 0 );
	}

	// update factored m1
	m1.ldlt_UpdateIncrement( w );

	m1.clearUpperTriangle();
	m2.clearUpperTriangle();

	if ( !m1.compare( m2, 1e-4f ) ) {
		Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::ldlt_UpdateIncrement failed" <<endl;
	}

	/*
		CMatXD::ldlt_UpdateDecrement
	*/

	for ( offset = 0; offset < size; offset++ ) {
		m1 = original;
		m2 = original;

		v.setSize( 6 );
		for ( int i = 0; i < size; i++ ) {
			v[i] = original[i][offset];
		}

		// factor m1
		m1.ldlt_Factor();

		// modify and factor m2
		m2.update_Decrement( offset );
		if ( !m2.ldlt_Factor() ) {
			SMF_ASSERT( 0 );
		}

		// update factors of m1
		m1.ldlt_UpdateDecrement( v, offset );

		if ( !m1.compare( m2, 1e-3f ) ) {
			Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::ldlt_UpdateDecrement failed" <<endl;
		}
	}

	/*
		CMatXD::ldlt_Inverse
	*/

	m2 = original;

	m2.ldlt_Factor();
	m2.ldlt_Inverse( m1 );
	m1 *= original;

	if ( !m1.isIdentity( 1e-4f ) ) {
		Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::ldlt_Inverse failed" <<endl;
	}

	/*
		CMatXD::eigen_solveSymmetricTriDiagonal
	*/

	m3 = original;
	m3.triDiagonal_ClearTriangles();
	m1 = m3;

	v.setSize( size );

	m1.eigen_solveSymmetricTriDiagonal( v );

	m3.transposeMultiply( m2, m1 );

	for ( int i = 0; i < size; i++ ) {
		for ( int j = 0; j < size; j++ ) {
			m1[i][j] *= v[j];
		}
	}

	if ( !m1.compare( m2, 1e-4f ) ) {
		Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::eigen_solveSymmetricTriDiagonal failed" <<endl;
	}

	/*
		CMatXD::eigen_solveSymmetric
	*/

	m3 = original;
	m1 = m3;

	v.setSize( size );

	m1.eigen_solveSymmetric( v );

	m3.transposeMultiply( m2, m1 );

	for ( int i = 0; i < size; i++ ) {
		for ( int j = 0; j < size; j++ ) {
			m1[i][j] *= v[j];
		}
	}

	if ( !m1.compare( m2, 1e-4f ) ) {
		Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::eigen_solveSymmetric failed" <<endl;
	}

	/*
		CMatXD::eigen_solve
	*/

	m3 = original;
	m1 = m3;

	v.setSize( size );
	w.setSize( size );

	m1.eigen_solve( v, w );

	m3.transposeMultiply( m2, m1 );

	for ( int i = 0; i < size; i++ ) {
		for ( int j = 0; j < size; j++ ) {
			m1[i][j] *= v[j];
		}
	}

	if ( !m1.compare( m2, 1e-4f ) ) {
		Debug::debug(Debug::error,__FUNCTION__)<<"CMatXD::eigen_solve failed" <<endl;
	}
}
} //END MATH
}//end SMF
