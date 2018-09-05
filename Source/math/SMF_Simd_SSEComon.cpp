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

#include "math/SMF_Math.h"
#include "math/SMF_Matriz.h"
#include "math/SMF_Vector.h"
#include "math/SMF_Simd.h"
#include "math/SMF_SimdGeneric.h"
#include "math/SMF_SimdMMX.h"
#include "math/SMF_SimdSSE.h"
#include "math/SMF_SimdSSE2.h"
#include "math/SMF_SimdSSE3.h"
#include "math/SMF_SimdSSE41.h"
#include "math/SMF_SimdAVX.h"
#include "math/SMF_JointTransform.h"
#include "util/SMF_Debug.h"
#include "geometry/all.h"

namespace SMF {
namespace MATH{
/// The returned SP FP contains x+y+z+w in all channels of the vector.
static  sf_m128 sum_xyzw_ps(sf_m128 m)
{
	m = _mm_hadd_ps(m, m); // m = (x+y, z+w, x+y, z+w).
	m = _mm_hadd_ps(m, m); // m = (x+y+z+w, x+y+z+w, x+y+z+w, x+y+z+w).
	return m; // Each index of the output will contain the sum x+y+z+w.
}
/// The dot product is stored in each channel of the returned vector.
static sf_m128 dot4_ps(sf_m128 a, sf_m128 b)
{

//#ifdef MATH_SSE41 // If we have SSE 4.1, we can use the dpps (dot product) instruction, _mm_dp_ps intrinsic.
//	return _mm_dp_ps(a, b, 0xFF); // Choose to multiply x, y, z and w (0xF0 = 1111 0000), and store the output to all indices (0x0F == 0000 1111).
//#else // Otherwise, use SSE3 haddps or SSE1 with individual shuffling.
	return sum_xyzw_ps(_mm_mul_ps(a, b));
//#endif
}


sf_m128 CSIMD_SSE::mat4x4_mul_sse(const sf_m128 *matrix, sf_m128 vector)
{
	sf_m128 x = _mm_mul_ps(matrix[0], vector);
	sf_m128 y = _mm_mul_ps(matrix[1], vector);
	sf_m128 z = _mm_mul_ps(matrix[2], vector);
	sf_m128 w = _mm_mul_ps(matrix[3], vector);
	_MM_TRANSPOSE4_PS(x, y, z, w); // Contains 2x unpacklo's, 2x unpackhi's, 2x movelh's and 2x movehl's. (or 8 shuffles, depending on the compiler)

	return _mm_add_ps(_mm_add_ps(x, y), _mm_add_ps(z, w));
}



sf_m128 CSIMD_SSE::colmajor_mat4x4_mul_sse(const sf_m128 *matrix, sf_m128 vector)
{
	sf_m128 x = _mm_shuffle_ps(vector, vector, _MM_SHUFFLE(0,0,0,0));
	sf_m128 y = _mm_shuffle_ps(vector, vector, _MM_SHUFFLE(1,1,1,1));
	sf_m128 z = _mm_shuffle_ps(vector, vector, _MM_SHUFFLE(2,2,2,2));
	sf_m128 w = _mm_shuffle_ps(vector, vector, _MM_SHUFFLE(3,3,3,3));
	x = _mm_mul_ps(x, matrix[0]);
	y = _mm_mul_ps(y, matrix[1]);
	z = _mm_mul_ps(z, matrix[2]);
	w = _mm_mul_ps(w, matrix[3]);

	return _mm_add_ps(_mm_add_ps(x, y), _mm_add_ps(z, w));
}


sf_m128 CSIMD_SSE::mat3x4_mul_sse(const sf_m128 *matrix, sf_m128 vector)
{
	sf_m128 x = dot4_ps(matrix[0], vector);
	sf_m128 y = dot4_ps(matrix[1], vector);
	sf_m128 z = dot4_ps(matrix[2], vector);

	// Take the 'w' component of the vector unmodified.
	sf_m128 xy = _mm_movelh_ps(x, y); // xy = [ _, y, _, x]
	sf_m128 zw = _mm_movehl_ps(vector, z); // zw = [ w, _, z, _]
	return _mm_shuffle_ps(xy, zw, _MM_SHUFFLE(3, 1, 2, 0)); // ret = [w, z, y, x]
}

CVec3D CSIMD_SSE::mat3x4_mul_vec(const sf_m128 *matrix, sf_m128 vector)
{
	sf_m128 x = dot4_ps(matrix[0], vector);
	sf_m128 y = dot4_ps(matrix[1], vector);
	sf_m128 z = dot4_ps(matrix[2], vector);

	return CVec3D(M128_TO_FLOAT(x), M128_TO_FLOAT(y), M128_TO_FLOAT(z));
}

void CSIMD_SSE::mat4x4_mul_sse(sf_m128 *out, const sf_m128 *m1, const sf_m128 *m2)
{
	sf_m128 s0 = shuffle1_ps(m1[0], _MM_SHUFFLE(0,0,0,0));
	sf_m128 s1 = shuffle1_ps(m1[0], _MM_SHUFFLE(1,1,1,1));
	sf_m128 s2 = shuffle1_ps(m1[0], _MM_SHUFFLE(2,2,2,2));
	sf_m128 s3 = shuffle1_ps(m1[0], _MM_SHUFFLE(3,3,3,3));
	sf_m128 r0 = _mm_mul_ps(s0, m2[0]);
	sf_m128 r1 = _mm_mul_ps(s1, m2[1]);
	sf_m128 r2 = _mm_mul_ps(s2, m2[2]);
	sf_m128 r3 = _mm_mul_ps(s3, m2[3]);
	out[0] = _mm_add_ps(_mm_add_ps(r0, r1), _mm_add_ps(r2, r3));

	s0 = shuffle1_ps(m1[1], _MM_SHUFFLE(0,0,0,0));
	s1 = shuffle1_ps(m1[1], _MM_SHUFFLE(1,1,1,1));
	s2 = shuffle1_ps(m1[1], _MM_SHUFFLE(2,2,2,2));
	s3 = shuffle1_ps(m1[1], _MM_SHUFFLE(3,3,3,3));
	r0 = _mm_mul_ps(s0, m2[0]);
	r1 = _mm_mul_ps(s1, m2[1]);
	r2 = _mm_mul_ps(s2, m2[2]);
	r3 = _mm_mul_ps(s3, m2[3]);
	out[1] = _mm_add_ps(_mm_add_ps(r0, r1), _mm_add_ps(r2, r3));

	s0 = shuffle1_ps(m1[2], _MM_SHUFFLE(0,0,0,0));
	s1 = shuffle1_ps(m1[2], _MM_SHUFFLE(1,1,1,1));
	s2 = shuffle1_ps(m1[2], _MM_SHUFFLE(2,2,2,2));
	s3 = shuffle1_ps(m1[2], _MM_SHUFFLE(3,3,3,3));
	r0 = _mm_mul_ps(s0, m2[0]);
	r1 = _mm_mul_ps(s1, m2[1]);
	r2 = _mm_mul_ps(s2, m2[2]);
	r3 = _mm_mul_ps(s3, m2[3]);
	out[2] = _mm_add_ps(_mm_add_ps(r0, r1), _mm_add_ps(r2, r3));

	s0 = shuffle1_ps(m1[3], _MM_SHUFFLE(0,0,0,0));
	s1 = shuffle1_ps(m1[3], _MM_SHUFFLE(1,1,1,1));
	s2 = shuffle1_ps(m1[3], _MM_SHUFFLE(2,2,2,2));
	s3 = shuffle1_ps(m1[3], _MM_SHUFFLE(3,3,3,3));
	r0 = _mm_mul_ps(s0, m2[0]);
	r1 = _mm_mul_ps(s1, m2[1]);
	r2 = _mm_mul_ps(s2, m2[2]);
	r3 = _mm_mul_ps(s3, m2[3]);
	out[3] = _mm_add_ps(_mm_add_ps(r0, r1), _mm_add_ps(r2, r3));
}

void CSIMD_SSE::mat3x4_mul_sse(sf_m128 *out, const sf_m128 *m1, const sf_m128 *m2)
{
	const sf_m128 m2_3 = MATH::SET_PS(1.f, 0.f, 0.f, 0.f);

	sf_m128 s0 = shuffle1_ps(m1[0], _MM_SHUFFLE(0,0,0,0));
	sf_m128 s1 = shuffle1_ps(m1[0], _MM_SHUFFLE(1,1,1,1));
	sf_m128 s2 = shuffle1_ps(m1[0], _MM_SHUFFLE(2,2,2,2));
	sf_m128 s3 = shuffle1_ps(m1[0], _MM_SHUFFLE(3,3,3,3));
	sf_m128 r0 = _mm_mul_ps(s0, m2[0]);
	sf_m128 r1 = _mm_mul_ps(s1, m2[1]);
	sf_m128 r2 = _mm_mul_ps(s2, m2[2]);
	sf_m128 r3 = _mm_mul_ps(s3, m2_3);
	out[0] = _mm_add_ps(_mm_add_ps(r0, r1), _mm_add_ps(r2, r3));

	s0 = shuffle1_ps(m1[1], _MM_SHUFFLE(0,0,0,0));
	s1 = shuffle1_ps(m1[1], _MM_SHUFFLE(1,1,1,1));
	s2 = shuffle1_ps(m1[1], _MM_SHUFFLE(2,2,2,2));
	s3 = shuffle1_ps(m1[1], _MM_SHUFFLE(3,3,3,3));
	r0 = _mm_mul_ps(s0, m2[0]);
	r1 = _mm_mul_ps(s1, m2[1]);
	r2 = _mm_mul_ps(s2, m2[2]);
	r3 = _mm_mul_ps(s3, m2_3);
	out[1] = _mm_add_ps(_mm_add_ps(r0, r1), _mm_add_ps(r2, r3));

	s0 = shuffle1_ps(m1[2], _MM_SHUFFLE(0,0,0,0));
	s1 = shuffle1_ps(m1[2], _MM_SHUFFLE(1,1,1,1));
	s2 = shuffle1_ps(m1[2], _MM_SHUFFLE(2,2,2,2));
	s3 = shuffle1_ps(m1[2], _MM_SHUFFLE(3,3,3,3));
	r0 = _mm_mul_ps(s0, m2[0]);
	r1 = _mm_mul_ps(s1, m2[1]);
	r2 = _mm_mul_ps(s2, m2[2]);
	r3 = _mm_mul_ps(s3, m2_3);
	out[2] = _mm_add_ps(_mm_add_ps(r0, r1), _mm_add_ps(r2, r3));
}

float CSIMD_SSE::mat4x4_inverse(const CMat4D *matrix, CMat4D *out)
{
	const sf_m128  *mat = matrix->Row(0).toM128Ptr();


	sf_m128 f1 = MAT_COFACTOR(mat, 3, 2);
	sf_m128 f2 = MAT_COFACTOR(mat, 3, 1);
	sf_m128 f3 = MAT_COFACTOR(mat, 2, 1);
	sf_m128 f4 = MAT_COFACTOR(mat, 3, 0);
	sf_m128 f5 = MAT_COFACTOR(mat, 2, 0);
	sf_m128 f6 = MAT_COFACTOR(mat, 1, 0);
	sf_m128 v1 = shuffle1_ps(_mm_shuffle_ps(mat[1], mat[0], _MM_SHUFFLE(0,0,0,0)), _MM_SHUFFLE(2,2,2,0));
	sf_m128 v2 = shuffle1_ps(_mm_shuffle_ps(mat[1], mat[0], _MM_SHUFFLE(1,1,1,1)), _MM_SHUFFLE(2,2,2,0));
	sf_m128 v3 = shuffle1_ps(_mm_shuffle_ps(mat[1], mat[0], _MM_SHUFFLE(2,2,2,2)), _MM_SHUFFLE(2,2,2,0));
	sf_m128 v4 = shuffle1_ps(_mm_shuffle_ps(mat[1], mat[0], _MM_SHUFFLE(3,3,3,3)), _MM_SHUFFLE(2,2,2,0));
	const sf_m128 s1 = _mm_set_ps(-0.0f,  0.0f, -0.0f,  0.0f);
	const sf_m128 s2 = _mm_set_ps( 0.0f, -0.0f,  0.0f, -0.0f);
	sf_m128 r1 = _mm_xor_ps(s1, _mm_add_ps(_mm_sub_ps(_mm_mul_ps(v2, f1), _mm_mul_ps(v3, f2)), _mm_mul_ps(v4, f3)));
	sf_m128 r2 = _mm_xor_ps(s2, _mm_add_ps(_mm_sub_ps(_mm_mul_ps(v1, f1), _mm_mul_ps(v3, f4)), _mm_mul_ps(v4, f5)));
	sf_m128 r3 = _mm_xor_ps(s1, _mm_add_ps(_mm_sub_ps(_mm_mul_ps(v1, f2), _mm_mul_ps(v2, f4)), _mm_mul_ps(v4, f6)));
	sf_m128 r4 = _mm_xor_ps(s2, _mm_add_ps(_mm_sub_ps(_mm_mul_ps(v1, f3), _mm_mul_ps(v2, f5)), _mm_mul_ps(v3, f6)));
	sf_m128 det = dot4_ps(mat[0], _mm_movelh_ps(_mm_unpacklo_ps(r1, r2), _mm_unpacklo_ps(r3, r4)));
	sf_m128 rcp = _mm_rcp_ps(det);

	out->Row(0).v = _mm_mul_ps(r1, rcp);
	out->Row(1).v = _mm_mul_ps(r2, rcp);
	out->Row(2).v = _mm_mul_ps(r3, rcp);
	out->Row(3).v = _mm_mul_ps(r4, rcp);
	return M128_TO_FLOAT(det);
}
static sf_m128 colmajor_mat4x4_muldir_sse(const sf_m128 *matrix, sf_m128 vector)
{
	sf_m128 x = _mm_shuffle_ps(vector, vector, _MM_SHUFFLE(0,0,0,0));
	sf_m128 y = _mm_shuffle_ps(vector, vector, _MM_SHUFFLE(1,1,1,1));
	sf_m128 z = _mm_shuffle_ps(vector, vector, _MM_SHUFFLE(2,2,2,2));
	x = _mm_mul_ps(x, matrix[0]);
	y = _mm_mul_ps(y, matrix[1]);
	z = _mm_mul_ps(z, matrix[2]);

	return _mm_add_ps(_mm_add_ps(x, y), z);
}

/// Inverts a 3x4 affine transformation matrix (in row-major format) that only consists of rotation (+possibly mirroring) and translation.
void CSIMD_SSE::mat3x4_inverse_orthonormal(sf_m128 *mat, sf_m128 *out)
{
	// mat[0]: [tx,02,01,00]
	// mat[1]: [ty,12,11,10]
	// mat[2]: [tz,22,21,20]
	// mat[3]: assumed to be [1,0,0,0] - not read.

	// First get the translation part (tx,ty,tz) from the original matrix,
	// and compute T=-M^(-1).
	sf_m128 tmp1 = _mm_unpackhi_ps(mat[0], mat[1]);                      // [ty,tx,12,02]
	sf_m128 xyz = _mm_shuffle_ps(tmp1, mat[2], _MM_SHUFFLE(3, 3, 3, 2)); // [ _,tz,ty,tx]
	sf_m128 vec = negate_ps(colmajor_mat4x4_muldir_sse(mat, xyz));      // [ _,Tz,Ty,Tx]

	sf_m128 tmp0 = _mm_unpacklo_ps(mat[0], mat[1]); // [11,01,10,00]
	//     tmp1 = computed already above              [ty,tx,12,02]
	sf_m128 tmp2 = _mm_unpacklo_ps(mat[2], vec);    // [Ty,21,Tx,20]
	sf_m128 tmp3 = _mm_unpackhi_ps(mat[2], vec);    // [ _,23,Tz,22]

	out[0] = _mm_movelh_ps(tmp0, tmp2);            // [Tx,20,10,00]
	out[1] = _mm_movehl_ps(tmp2, tmp0);            // [Ty,21,11,01]
	out[2] = _mm_movelh_ps(tmp1, tmp3);            // [Tz,22 12,02]
 // out[3] = assumed to be [1,0,0,0] - no need to write back.
}


sf_m128 CSIMD_SSE:: newtonRhapsonRecipStep(sf_m128 recip, sf_m128 estimate)
{
	// Do one iteration of Newton-Rhapson:
	// e_n = 2*e - x*e^2
	sf_m128 e2 = _mm_mul_ps(estimate, estimate);
	return _mm_sub_ps(_mm_add_ps(estimate, estimate), _mm_mul_ps(recip, e2));
}

sf_m128 CSIMD_SSE:: newtonRhapsonRecip(sf_m128 recip)
{
	sf_m128 estimate = _mm_rcp_ps(recip);
	return  newtonRhapsonRecipStep(recip, estimate);
}

/// Computes the determinant of a 4x4 matrix.
float CSIMD_SSE::mat4x4_determinant(const CMat4D *matrix)
{
	const sf_m128 *row = reinterpret_cast<const sf_m128 *>(&matrix);
	sf_m128 s = shuffle1_ps( newtonRhapsonRecip(row[0]), _MM_SHUFFLE(0,0,0,0));
	// row[0].x has a factor of the final determinant.
	sf_m128 row0 = _mm_mul_ps(s, row[0]);
	s = shuffle1_ps(row[1], _MM_SHUFFLE(0,0,0,0));
	sf_m128 row1 = _mm_sub_ps(row[1], _mm_mul_ps(s, row0));
	s = shuffle1_ps(row[2], _MM_SHUFFLE(0,0,0,0));
	sf_m128 row2 = _mm_sub_ps(row[2], _mm_mul_ps(s, row0));
	s = shuffle1_ps(row[3], _MM_SHUFFLE(0,0,0,0));
	sf_m128 row3 = _mm_sub_ps(row[3], _mm_mul_ps(s, row0));

	// row1.y has a factor of the final determinant.
	s = shuffle1_ps( newtonRhapsonRecip(row1), _MM_SHUFFLE(1,1,1,1));
	sf_m128 row1_1 = _mm_mul_ps(s, row1);
	s = shuffle1_ps(row2, _MM_SHUFFLE(1,1,1,1));
	sf_m128 row2_1 = _mm_sub_ps(row2, _mm_mul_ps(s, row1_1));
	s = shuffle1_ps(row3, _MM_SHUFFLE(1,1,1,1));
	sf_m128 row3_1 = _mm_sub_ps(row3, _mm_mul_ps(s, row1_1));

	// Now we are left with a 2x2 matrix in row2_1.zw and row3_1.zw.
	// D = row2_1.z * row3_1.w - row2_1.w * row3_1.z.
	sf_m128 r1 = shuffle1_ps(row2_1, _MM_SHUFFLE(2,3,1,0));
	sf_m128 r = _mm_mul_ps(r1, row3_1);
	sf_m128 a = shuffle1_ps(r, _MM_SHUFFLE(3,3,3,3));
	sf_m128 b = shuffle1_ps(r, _MM_SHUFFLE(2,2,2,2));
	sf_m128 d1 = _mm_sub_ss(a, b);
	sf_m128 d2 = row[0];
	sf_m128 d3 = shuffle1_ps(row1, _MM_SHUFFLE(1,1,1,1));
	sf_m128 d = _mm_mul_ss(d1, _mm_mul_ss(d2, d3));
	return M128_TO_FLOAT(d);
}

/// Computes the determinant of a 3x4 matrix stored in row-major format. (Treated as a square matrix with last row [0,0,0,1])
float CSIMD_SSE::mat3x4_determinant(const sf_m128 *row)
{
	sf_m128 s = shuffle1_ps( newtonRhapsonRecip(row[0]), _MM_SHUFFLE(0,0,0,0));
	// row[0].x has a factor of the final determinant.
	sf_m128 row0 = _mm_mul_ps(s, row[0]);
	s = shuffle1_ps(row[1], _MM_SHUFFLE(0,0,0,0));
	sf_m128 row1 = _mm_sub_ps(row[1], _mm_mul_ps(s, row0));
	s = shuffle1_ps(row[2], _MM_SHUFFLE(0,0,0,0));
	sf_m128 row2 = _mm_sub_ps(row[2], _mm_mul_ps(s, row0));

	// Now we are left with a 2x2 matrix in row1.yz and row2.yz.
	// D = row1.y * row2.z - row2.y * row1.z.
	sf_m128 r1 = shuffle1_ps(row1, _MM_SHUFFLE(3,1,2,0));
	sf_m128 r = _mm_mul_ps(r1, row2);
	sf_m128 a = shuffle1_ps(r, _MM_SHUFFLE(2,2,2,2));
	sf_m128 b = shuffle1_ps(r, _MM_SHUFFLE(1,1,1,1));
	sf_m128 d1 = _mm_sub_ss(a, b);
	sf_m128 d2 = row[0];
	sf_m128 d = _mm_mul_ss(d1, d2);
	return M128_TO_FLOAT(d);
}

void CSIMD_SSE::mat3x4_transpose(const sf_m128 *src, sf_m128 *dst)
{
	sf_m128 src3 = _mm_setzero_ps(); // w component should be 1, but since it won't get stored, it doesn't matter, so we can just create zeros.
	sf_m128 tmp0 = _mm_unpacklo_ps(src[0], src[1]);
	sf_m128 tmp2 = _mm_unpacklo_ps(src[2], src3);
	sf_m128 tmp1 = _mm_unpackhi_ps(src[0], src[1]);
	sf_m128 tmp3 = _mm_unpackhi_ps(src[2], src3);
	dst[0] = _mm_movelh_ps(tmp0, tmp2);
	dst[1] = _mm_movehl_ps(tmp2, tmp0);
	dst[2] = _mm_movelh_ps(tmp1, tmp3);
}

void CSIMD_SSE::mat4x4_mul_dpps(sf_m128 *out, const sf_m128 *m1, const sf_m128 *m2)
{
	// Transpose m2:
	// m2[0] = [ 03, 02, 01, 00 ]     [ 30, 20, 10, 00 ]
	// m2[1] = [ 13, 12, 11, 10 ] --> [ 31, 21, 11, 01 ]
	// m2[2] = [ 23, 22, 21, 20 ] --> [ 32, 22, 12, 02 ]
	//         [ 33, 32, 31, 30 ]     [ 33, 23, 13, 03 ]

	sf_m128 low1 = _mm_movelh_ps(m2[0], m2[1]); // = [ 11, 10, 01, 00 ]
	sf_m128 low2 = _mm_movelh_ps(m2[2], m2[3]); // = [ 31, 30, 21, 20 ]
	sf_m128 hi1 = _mm_movehl_ps(m2[1], m2[0]);  // = [ 13, 12, 03, 02 ]
	sf_m128 hi2 = _mm_movehl_ps(m2[3], m2[2]);  // = [ 33, 32, 23, 22 ]

	sf_m128 row1 = _mm_shuffle_ps(low1, low2, _MM_SHUFFLE(2, 0, 2, 0)); // = [30, 20, 10, 00]
	sf_m128 row2 = _mm_shuffle_ps(low1, low2, _MM_SHUFFLE(3, 1, 3, 1)); // = [31, 21, 11, 01]
	sf_m128 row3 = _mm_shuffle_ps(hi1, hi2, _MM_SHUFFLE(2, 0, 2, 0));   // = [32, 22, 12, 02]
	sf_m128 row4 = _mm_shuffle_ps(hi1, hi2, _MM_SHUFFLE(3, 1, 3, 1));   // = [33, 23, 13, 03]

	sf_m128 _00 = dot4_ps(m1[0], row1);
	sf_m128 _01 = dot4_ps(m1[0], row2);
	sf_m128 _02 = dot4_ps(m1[0], row3);
	sf_m128 _03 = dot4_ps(m1[0], row4);
	out[0] = pack_4ss_to_ps(_00, _01, _02, _03);

	sf_m128 _10 = dot4_ps(m1[1], row1);
	sf_m128 _11 = dot4_ps(m1[1], row2);
	sf_m128 _12 = dot4_ps(m1[1], row3);
	sf_m128 _13 = dot4_ps(m1[1], row4);
	out[1] = pack_4ss_to_ps(_10, _11, _12, _13);

	sf_m128 _20 = dot4_ps(m1[2], row1);
	sf_m128 _21 = dot4_ps(m1[2], row2);
	sf_m128 _22 = dot4_ps(m1[2], row3);
	sf_m128 _23 = dot4_ps(m1[2], row4);
	out[2] = pack_4ss_to_ps(_20, _21, _22, _23);

	sf_m128 _30 = dot4_ps(m1[3], row1);
	sf_m128 _31 = dot4_ps(m1[3], row2);
	sf_m128 _32 = dot4_ps(m1[3], row3);
	sf_m128 _33 = dot4_ps(m1[3], row4);
	out[3] = pack_4ss_to_ps(_30, _31, _32, _33);
}

void CSIMD_SSE::mat4x4_mul_dpps_2(sf_m128 *out, const sf_m128 *m1, const sf_m128 *m2)
{
	// Transpose m2:
	// m2[0] = [ 03, 02, 01, 00 ]     [ 30, 20, 10, 00 ]
	// m2[1] = [ 13, 12, 11, 10 ] --> [ 31, 21, 11, 01 ]
	// m2[2] = [ 23, 22, 21, 20 ] --> [ 32, 22, 12, 02 ]
	//         [ 33, 32, 31, 30 ]     [ 33, 23, 13, 03 ]
	sf_m128 row1 = m2[0];
	sf_m128 row2 = m2[1];
	sf_m128 row3 = m2[2];
	sf_m128 row4 = m2[3];
	_mm_transpose_matrix_intel(row1, row2, row3, row4);

	sf_m128 _00 = dot4_ps(m1[0], row1);
	sf_m128 _01 = dot4_ps(m1[0], row2);
	sf_m128 _02 = dot4_ps(m1[0], row3);
	sf_m128 _03 = dot4_ps(m1[0], row4);
	out[0] = pack_4ss_to_ps(_00, _01, _02, _03);

	sf_m128 _10 = dot4_ps(m1[1], row1);
	sf_m128 _11 = dot4_ps(m1[1], row2);
	sf_m128 _12 = dot4_ps(m1[1], row3);
	sf_m128 _13 = dot4_ps(m1[1], row4);
	out[1] = pack_4ss_to_ps(_10, _11, _12, _13);

	sf_m128 _20 = dot4_ps(m1[2], row1);
	sf_m128 _21 = dot4_ps(m1[2], row2);
	sf_m128 _22 = dot4_ps(m1[2], row3);
	sf_m128 _23 = dot4_ps(m1[2], row4);
	out[2] = pack_4ss_to_ps(_20, _21, _22, _23);

	sf_m128 _30 = dot4_ps(m1[3], row1);
	sf_m128 _31 = dot4_ps(m1[3], row2);
	sf_m128 _32 = dot4_ps(m1[3], row3);
	sf_m128 _33 = dot4_ps(m1[3], row4);
	out[3] = pack_4ss_to_ps(_30, _31, _32, _33);
}

void CSIMD_SSE::mat4x4_mul_dpps_3(sf_m128 *out, const sf_m128 *m1, const sf_m128 *m2)
{
	// Transpose m2:
	// m2[0] = [ 03, 02, 01, 00 ]     [ 30, 20, 10, 00 ]
	// m2[1] = [ 13, 12, 11, 10 ] --> [ 31, 21, 11, 01 ]
	// m2[2] = [ 23, 22, 21, 20 ] --> [ 32, 22, 12, 02 ]
	//         [ 33, 32, 31, 30 ]     [ 33, 23, 13, 03 ]

	sf_m128 low1 = _mm_movelh_ps(m2[0], m2[1]); // = [ 11, 10, 01, 00 ]
	sf_m128 low2 = _mm_movelh_ps(m2[2], m2[3]); // = [ 31, 30, 21, 20 ]
	sf_m128 hi1 = _mm_movehl_ps(m2[1], m2[0]);  // = [ 13, 12, 03, 02 ]
	sf_m128 hi2 = _mm_movehl_ps(m2[3], m2[2]);  // = [ 33, 32, 23, 22 ]

	sf_m128 row1 = _mm_shuffle_ps(low1, low2, _MM_SHUFFLE(2, 0, 2, 0)); // = [30, 20, 10, 00]
	sf_m128 row2 = _mm_shuffle_ps(low1, low2, _MM_SHUFFLE(3, 1, 3, 1)); // = [31, 21, 11, 01]
	sf_m128 row3 = _mm_shuffle_ps(hi1, hi2, _MM_SHUFFLE(2, 0, 2, 0));   // = [32, 22, 12, 02]
	sf_m128 row4 = _mm_shuffle_ps(hi1, hi2, _MM_SHUFFLE(3, 1, 3, 1));   // = [33, 23, 13, 03]

	sf_m128 _00 = dot4_ps(m1[0], row1);
	sf_m128 _01 = dot4_ps(m1[0], row2);
	sf_m128 _02 = dot4_ps(m1[0], row3);
	sf_m128 _03 = dot4_ps(m1[0], row4);

	sf_m128 xy = _mm_movelh_ps(_00, _01); // xy = [ _, y, _, x]
	sf_m128 zw = _mm_movelh_ps(_02, _03); // zw = [ _, w, _, z]
	out[0] = _mm_shuffle_ps(xy, zw, _MM_SHUFFLE(2, 0, 2, 0)); // ret = [w, z, y, x]

//	out[0] = pack_4ss_to_ps(_00, _01, _02, _03);

	sf_m128 _10 = dot4_ps(m1[1], row1);
	sf_m128 _11 = dot4_ps(m1[1], row2);
	sf_m128 _12 = dot4_ps(m1[1], row3);
	sf_m128 _13 = dot4_ps(m1[1], row4);

	sf_m128 xy2 = _mm_movelh_ps(_10, _11); // xy = [ _, y, _, x]
	sf_m128 zw2 = _mm_movelh_ps(_12, _13); // zw = [ _, w, _, z]
	out[1] = _mm_shuffle_ps(xy2, zw2, _MM_SHUFFLE(2, 0, 2, 0)); // ret = [w, z, y, x]

//	out[1] = pack_4ss_to_ps(_10, _11, _12, _13);

	sf_m128 _20 = dot4_ps(m1[2], row1);
	sf_m128 _21 = dot4_ps(m1[2], row2);
	sf_m128 _22 = dot4_ps(m1[2], row3);
	sf_m128 _23 = dot4_ps(m1[2], row4);

	sf_m128 xy3 = _mm_movelh_ps(_20, _21); // xy = [ _, y, _, x]
	sf_m128 zw3 = _mm_movelh_ps(_22, _23); // zw = [ _, w, _, z]
	out[2] = _mm_shuffle_ps(xy3, zw3, _MM_SHUFFLE(2, 0, 2, 0)); // ret = [w, z, y, x]

//	out[2] = pack_4ss_to_ps(_20, _21, _22, _23);

	sf_m128 _30 = dot4_ps(m1[3], row1);
	sf_m128 _31 = dot4_ps(m1[3], row2);
	sf_m128 _32 = dot4_ps(m1[3], row3);
	sf_m128 _33 = dot4_ps(m1[3], row4);

	sf_m128 xy4 = _mm_movelh_ps(_30, _31); // xy = [ _, y, _, x]
	sf_m128 zw4 = _mm_movelh_ps(_32, _33); // zw = [ _, w, _, z]
	out[3] = _mm_shuffle_ps(xy4, zw4, _MM_SHUFFLE(2, 0, 2, 0)); // ret = [w, z, y, x]

//	out[3] = pack_4ss_to_ps(_30, _31, _32, _33);
}

void CSIMD_SSE::quat_to_mat3x4(sf_m128 q, sf_m128 t, sf_m128 *m)
{
	// Constants:
	const sf_u32 sign = 0x80000000UL;
	const sf_m128 sseX0 = set_ps_hex(sign, sign, sign, 0);
	const sf_m128 sseX1 = set_ps_hex(sign, sign, 0, sign);
	sf_m128 one = _mm_set_ps(0, 0, 0, 1);

#if 0 // The original code converted a quaternion into an hybrid of rotation/translation (bug?)
	sf_m128 q2 = _mm_add_ps(q, q);                                 // [2w 2z 2y 2x]
	sf_m128 yxxy = shuffle1_ps(q, _MM_SHUFFLE(1, 0, 0, 1));        // [ y  x  x  y]
	sf_m128 yyzz2 = shuffle1_ps(q2, _MM_SHUFFLE(2, 2, 1, 1));      // [2z 2z 2y 2y]
	sf_m128 yy_xy_xz_yz_2 = _mm_mul_ps(yxxy, yyzz2);               // [2yz 2xz 2xy 2yy]

	sf_m128 zwww = shuffle1_ps(q, _MM_SHUFFLE(3, 3, 3, 2));        // [w w w z]
	sf_m128 zzyx2 = shuffle1_ps(q2, _MM_SHUFFLE(0, 1, 2, 2));      // [2x 2y 2z 2z]
	sf_m128 zz_wz_wy_wx_2 = _mm_mul_ps(zwww, zzyx2);               // [2xw 2yw 2zw 2zz]

	sf_m128 xx2 = _mm_mul_ss(q, q2);                               // [2xx]

	// Calculate last two elements of the third row.
	sf_m128 one_m_xx2 = _mm_sub_ss(one, xx2);                      // [0 0 0 1-2xx]
	sf_m128 one_m_xx_yy_2 = _mm_sub_ss(one_m_xx2, yy_xy_xz_yz_2);  // [0 0 0 1-2xx-2yy]
	sf_m128 one_m_xx_yy_2_0_tz_tw = _mm_shuffle_ps(one_m_xx_yy_2, t, _MM_SHUFFLE(3, 2, 1, 0)); // [tw tz 0 1-2xx-2yy]

	// Calculate first row
	sf_m128 m_yy_xy_xz_yz_2 = _mm_xor_ps(yy_xy_xz_yz_2, sseX0);     // [-2yz -2xz -2xy   2yy]
	sf_m128 m_zz_wz_wy_wx_2 = _mm_xor_ps(zz_wz_wy_wx_2, sseX1);     // [-2xw -2yw  2zw  -2zz]
	sf_m128 m_zz_one_wz_wy_wx_2 = _mm_add_ss(m_zz_wz_wy_wx_2, one); // [-2xw -2yw  2zw 1-2zz]
	sf_m128 first_row = _mm_sub_ps(m_zz_one_wz_wy_wx_2, m_yy_xy_xz_yz_2); // [2yz-2xw 2xz-2yw 2xy+2zw 1-2zz-2yy]
	m[0] = first_row;
	_mm_store_ss((float*)m+3, t);

	// Calculate second row
	sf_m128 s1 = _mm_move_ss(m_yy_xy_xz_yz_2, xx2);                // [-2yz -2xz -2xy 2xx]
	sf_m128 s2 = _mm_xor_ps(m_zz_one_wz_wy_wx_2, sseX0);           // [2xw 2yw -2zw 1-2zz]
	sf_m128 s3 = _mm_sub_ps(s2, s1);                               // [2xw+2yz 2yw+2xz 2xy-2zw 1-2zz-2xx]
	sf_m128 t_yzwx = shuffle1_ps(t, _MM_SHUFFLE(0, 3, 2, 1));      // [tx tw tz ty]
	sf_m128 second_row = shuffle1_ps(s3, _MM_SHUFFLE(2, 3, 0, 1)); // [2yw+2xz 2xw+2yz 1-2zz-2xx 2xy-2zw]
	m[1] = second_row;
	_mm_store_ss((float*)m+7, t_yzwx);

	// Calculate third row
	sf_m128 t1 = _mm_movehl_ps(first_row, second_row);             // [2yz-2xw 2xz-2yw 2yw+2xz 2xw+2yz]
	sf_m128 t2 = _mm_shuffle_ps(t1, one_m_xx_yy_2_0_tz_tw, _MM_SHUFFLE(2, 0, 3, 1)); // [tz 1-2xx-2yy 2yz-2xw 2yw+2xz]
	m[2] = t2;
#else
	sf_m128 q2 = _mm_add_ps(q, q);                                 // [2w 2z 2y 2x]
	sf_m128 yxxy = shuffle1_ps(q, _MM_SHUFFLE(1, 0, 0, 1));        // [ y  x  x  y]
	sf_m128 yyzz2 = shuffle1_ps(q2, _MM_SHUFFLE(2, 2, 1, 1));      // [2z 2z 2y 2y]
	sf_m128 yy_xy_xz_yz_2 = _mm_mul_ps(yxxy, yyzz2);               // [2yz 2xz 2xy 2yy]

	sf_m128 zwww = shuffle1_ps(q, _MM_SHUFFLE(3, 3, 3, 2));        // [w w w z]
	sf_m128 zzyx2 = shuffle1_ps(q2, _MM_SHUFFLE(0, 1, 2, 2));      // [2x 2y 2z 2z]
	sf_m128 zz_wz_wy_wx_2 = _mm_mul_ps(zwww, zzyx2);               // [2xw 2yw 2zw 2zz]

	sf_m128 xx2 = _mm_mul_ss(q, q2);                               // [2xx]

	// Calculate last two elements of the third row.
	sf_m128 one_m_xx2 = _mm_sub_ss(one, xx2);                      // [0 0 0 1-2xx]
	sf_m128 one_m_xx_yy_2 = _mm_sub_ss(one_m_xx2, yy_xy_xz_yz_2);  // [0 0 0 1-2xx-2yy]
	sf_m128 one_m_xx_yy_2_0_tz_tw = one_m_xx_yy_2;//_mm_shuffle_ps(one_m_xx_yy_2, t, _MM_SHUFFLE(3, 2, 1, 0)); // [tw tz 0 1-2xx-2yy]

	// Calculate first row
	sf_m128 m_yy_xy_xz_yz_2 = _mm_xor_ps(yy_xy_xz_yz_2, sseX0);     // [-2yz -2xz -2xy   2yy]
	sf_m128 m_zz_wz_wy_wx_2 = _mm_xor_ps(zz_wz_wy_wx_2, sseX1);     // [-2xw -2yw  2zw  -2zz]
	sf_m128 m_zz_one_wz_wy_wx_2 = _mm_add_ss(m_zz_wz_wy_wx_2, one); // [-2xw -2yw  2zw 1-2zz]
	sf_m128 first_row = _mm_sub_ps(m_zz_one_wz_wy_wx_2, m_yy_xy_xz_yz_2); // [2yz-2xw 2xz-2yw 2xy+2zw 1-2zz-2yy]

	// Calculate second row
	sf_m128 s1 = _mm_move_ss(m_yy_xy_xz_yz_2, xx2);                // [-2yz -2xz -2xy 2xx]
	sf_m128 s2 = _mm_xor_ps(m_zz_one_wz_wy_wx_2, sseX0);           // [2xw 2yw -2zw 1-2zz]
	sf_m128 s3 = _mm_sub_ps(s2, s1);                               // [2xw+2yz 2yw+2xz 2xy-2zw 1-2zz-2xx]
	sf_m128 second_row = shuffle1_ps(s3, _MM_SHUFFLE(2, 3, 0, 1)); // [2yw+2xz 2xw+2yz 1-2zz-2xx 2xy-2zw]

	// Calculate third row
	sf_m128 t1 = _mm_movehl_ps(first_row, second_row);             // [2yz-2xw 2xz-2yw 2yw+2xz 2xw+2yz]
	sf_m128 third_row = _mm_shuffle_ps(t1, one_m_xx_yy_2_0_tz_tw, _MM_SHUFFLE(2, 0, 3, 1)); // [0 1-2xx-2yy 2yz-2xw 2yw+2xz]

	sf_m128 tmp0 = _mm_unpacklo_ps(first_row, second_row);
	sf_m128 tmp2 = _mm_unpacklo_ps(third_row, t);
	sf_m128 tmp1 = _mm_unpackhi_ps(first_row, second_row);
	sf_m128 tmp3 = _mm_unpackhi_ps(third_row, t);
	m[0] = _mm_movelh_ps(tmp0, tmp2);
	m[1] = _mm_movehl_ps(tmp2, tmp0);
	m[2] = _mm_movelh_ps(tmp1, tmp3);
#endif
}


void CSIMD_SSE::quat_to_mat4x4(sf_m128 q, sf_m128 t, sf_m128 *m)
{
	quat_to_mat3x4(q, t, m);
	m[3] = _mm_set_ps(1.f, 0.f, 0.f, 0.f);
}

sf_m128 CSIMD_SSE::quat_transform_vec4(sf_m128 quat, sf_m128 vec)
{
	sf_m128 W = shuffle1_ps(quat, _MM_SHUFFLE(3,3,3,3));

//	sf_m128 qxv = cross_ps(q, vec.v);
	sf_m128 a_xzy = shuffle1_ps(quat, _MM_SHUFFLE(3, 0, 2, 1)); // a_xzy = [a.w, a.x, a.z, a.y]
	sf_m128 b_yxz = shuffle1_ps(vec, _MM_SHUFFLE(3, 1, 0, 2)); // b_yxz = [b.w, b.y, b.x, b.z]
	sf_m128 a_yxz = shuffle1_ps(quat, _MM_SHUFFLE(3, 1, 0, 2)); // a_yxz = [a.w, a.y, a.x, a.z]
	sf_m128 b_xzy = shuffle1_ps(vec, _MM_SHUFFLE(3, 0, 2, 1)); // b_xzy = [b.w, b.x, b.z, b.y]
	sf_m128 x = _mm_mul_ps(a_xzy, b_yxz); // [a.w*b.w, a.x*b.y, a.z*b.x, a.y*b.z]
	sf_m128 y = _mm_mul_ps(a_yxz, b_xzy); // [a.w*b.w, a.y*b.x, a.x*b.z, a.z*b.y]
	sf_m128 qxv = _mm_sub_ps(x, y); // [0, a.x*b.y - a.y*b.x, a.z*b.x - a.x*b.z, a.y*b.z - a.z*b.y]

	sf_m128 Wv = _mm_mul_ps(W, vec);
	sf_m128 s = _mm_add_ps(qxv, Wv);

//	s = cross_ps(q, s);
	sf_m128 s_yxz = shuffle1_ps(s, _MM_SHUFFLE(3, 1, 0, 2)); // b_yxz = [b.w, b.y, b.x, b.z]
	sf_m128 s_xzy = shuffle1_ps(s, _MM_SHUFFLE(3, 0, 2, 1)); // b_xzy = [b.w, b.x, b.z, b.y]
	x = _mm_mul_ps(a_xzy, s_yxz); // [a.w*b.w, a.x*b.y, a.z*b.x, a.y*b.z]
	y = _mm_mul_ps(a_yxz, s_xzy); // [a.w*b.w, a.y*b.x, a.x*b.z, a.z*b.y]
	s = _mm_sub_ps(x, y); // [0, a.x*b.y - a.y*b.x, a.z*b.x - a.x*b.z, a.y*b.z - a.z*b.y]

	s = _mm_add_ps(s, s);
	s = _mm_add_ps(s, vec);
	return s;
}


sf_m128 CSIMD_SSE::quat_mul_quat(sf_m128 q1, sf_m128 q2)
{
/*	return Quat(x*r.w + y*r.z - z*r.y + w*r.x,
	           -x*r.z + y*r.w + z*r.x + w*r.y,
	            x*r.y - y*r.x + z*r.w + w*r.z,
	           -x*r.x - y*r.y - z*r.z + w*r.w); */

	const sf_m128 signx = set_ps_hex(0x80000000u, 0, 0x80000000u, 0); // [- + - +]
	const sf_m128 signy = shuffle1_ps(signx, _MM_SHUFFLE(3,3,0,0));   // [- - + +]
	const sf_m128 signz = shuffle1_ps(signx, _MM_SHUFFLE(3,0,0,3));   // [- + + -]

	sf_m128 X = _mm_xor_ps(signx, shuffle1_ps(q1, _MM_SHUFFLE(0,0,0,0)));
	sf_m128 Y = _mm_xor_ps(signy, shuffle1_ps(q1, _MM_SHUFFLE(1,1,1,1)));
	sf_m128 Z = _mm_xor_ps(signz, shuffle1_ps(q1, _MM_SHUFFLE(2,2,2,2)));
	sf_m128 W = shuffle1_ps(q1, _MM_SHUFFLE(3,3,3,3));

	sf_m128 r1 = shuffle1_ps(q2, _MM_SHUFFLE(0, 1, 2, 3)); // [x,y,z,w]
	sf_m128 r2 = shuffle1_ps(q2, _MM_SHUFFLE(1, 0, 3, 2)); // [y,x,w,z]
	sf_m128 r3 = shuffle1_ps(q2, _MM_SHUFFLE(2, 3, 0, 1)); // [z,w,x,y]
	// sf_m128 r4 = q2;

	return _mm_add_ps(_mm_add_ps(_mm_mul_ps(X, r1), _mm_mul_ps(Y, r2)),
	                  _mm_add_ps(_mm_mul_ps(Z, r3), _mm_mul_ps(W, q2)));

}


sf_m128 CSIMD_SSE::quat_div_quat(sf_m128 q1, sf_m128 q2)
{
/*	return Quat(x*r.w - y*r.z + z*r.y - w*r.x,
	            x*r.z + y*r.w - z*r.x - w*r.y,
	           -x*r.y + y*r.x + z*r.w - w*r.z,
	            x*r.x + y*r.y + z*r.z + w*r.w); */

	const sf_m128 signx = set_ps_hex(0x80000000u, 0, 0x80000000u, 0); // [- + - +]
	const sf_m128 signy = shuffle1_ps(signx, _MM_SHUFFLE(3,3,0,0));   // [- - + +]
	const sf_m128 signz = shuffle1_ps(signx, _MM_SHUFFLE(3,0,0,3));   // [- + + -]

	sf_m128 X = _mm_xor_ps(signx, shuffle1_ps(q1, _MM_SHUFFLE(0,0,0,0)));
	sf_m128 Y = _mm_xor_ps(signy, shuffle1_ps(q1, _MM_SHUFFLE(1,1,1,1)));
	sf_m128 Z = _mm_xor_ps(signz, shuffle1_ps(q1, _MM_SHUFFLE(2,2,2,2)));
	sf_m128 W = shuffle1_ps(q1, _MM_SHUFFLE(3,3,3,3));

	q2 = negate3_ps(q2);
	sf_m128 r1 = shuffle1_ps(q2, _MM_SHUFFLE(0, 1, 2, 3)); // [x,y,z,w]
	sf_m128 r2 = shuffle1_ps(q2, _MM_SHUFFLE(1, 0, 3, 2)); // [y,x,w,z]
	sf_m128 r3 = shuffle1_ps(q2, _MM_SHUFFLE(2, 3, 0, 1)); // [z,w,x,y]
	// sf_m128 r4 = q2;

	return _mm_add_ps(_mm_add_ps(_mm_mul_ps(X, r1), _mm_mul_ps(Y, r2)),
	                  _mm_add_ps(_mm_mul_ps(Z, r3), _mm_mul_ps(W, q2)));
}


void CSIMD_SSE2::mat4x4_mul_sse(sf_m128 *out, const sf_m128 *m1, const sf_m128 *m2)
{
	sf_m128 s0 = _mm_shuffle_ps(m1[0], m1[0], _MM_SHUFFLE(0,0,0,0));
	sf_m128 s1 = _mm_shuffle_ps(m1[0], m1[0], _MM_SHUFFLE(1,1,1,1));
	sf_m128 s2 = _mm_shuffle_ps(m1[0], m1[0], _MM_SHUFFLE(2,2,2,2));
	sf_m128 s3 = _mm_shuffle_ps(m1[0], m1[0], _MM_SHUFFLE(3,3,3,3));
	sf_m128 r0 = _mm_mul_ps(s0, m2[0]);
	sf_m128 r1 = _mm_mul_ps(s1, m2[1]);
	sf_m128 r2 = _mm_mul_ps(s2, m2[2]);
	sf_m128 r3 = _mm_mul_ps(s3, m2[3]);
	out[0] = _mm_add_ps(_mm_add_ps(r0, r1), _mm_add_ps(r2, r3));

	s0 = _mm_shuffle_ps(m1[1], m1[1], _MM_SHUFFLE(0,0,0,0));
	s1 = _mm_shuffle_ps(m1[1], m1[1], _MM_SHUFFLE(1,1,1,1));
	s2 = _mm_shuffle_ps(m1[1], m1[1], _MM_SHUFFLE(2,2,2,2));
	s3 = _mm_shuffle_ps(m1[1], m1[1], _MM_SHUFFLE(3,3,3,3));
	r0 = _mm_mul_ps(s0, m2[0]);
	r1 = _mm_mul_ps(s1, m2[1]);
	r2 = _mm_mul_ps(s2, m2[2]);
	r3 = _mm_mul_ps(s3, m2[3]);
	out[1] = _mm_add_ps(_mm_add_ps(r0, r1), _mm_add_ps(r2, r3));

	s0 = _mm_shuffle_ps(m1[2], m1[2], _MM_SHUFFLE(0,0,0,0));
	s1 = _mm_shuffle_ps(m1[2], m1[2], _MM_SHUFFLE(1,1,1,1));
	s2 = _mm_shuffle_ps(m1[2], m1[2], _MM_SHUFFLE(2,2,2,2));
	s3 = _mm_shuffle_ps(m1[2], m1[2], _MM_SHUFFLE(3,3,3,3));
	r0 = _mm_mul_ps(s0, m2[0]);
	r1 = _mm_mul_ps(s1, m2[1]);
	r2 = _mm_mul_ps(s2, m2[2]);
	r3 = _mm_mul_ps(s3, m2[3]);
	out[2] = _mm_add_ps(_mm_add_ps(r0, r1), _mm_add_ps(r2, r3));

	s0 = _mm_shuffle_ps(m1[3], m1[3], _MM_SHUFFLE(0,0,0,0));
	s1 = _mm_shuffle_ps(m1[3], m1[3], _MM_SHUFFLE(1,1,1,1));
	s2 = _mm_shuffle_ps(m1[3], m1[3], _MM_SHUFFLE(2,2,2,2));
	s3 = _mm_shuffle_ps(m1[3], m1[3], _MM_SHUFFLE(3,3,3,3));
	r0 = _mm_mul_ps(s0, m2[0]);
	r1 = _mm_mul_ps(s1, m2[1]);
	r2 = _mm_mul_ps(s2, m2[2]);
	r3 = _mm_mul_ps(s3, m2[3]);
	out[3] = _mm_add_ps(_mm_add_ps(r0, r1), _mm_add_ps(r2, r3));
}


 sf_m128 CSIMD_SSE2::colmajor_mat4x4_mul_sse(const sf_m128 *matrix, sf_m128 vector)
{
	sf_m128 x = shuffle2_ps(vector, _MM_SHUFFLE(0,0,0,0));
	sf_m128 y = shuffle2_ps(vector, _MM_SHUFFLE(1,1,1,1));
	sf_m128 z = shuffle2_ps(vector, _MM_SHUFFLE(2,2,2,2));
	sf_m128 w = shuffle2_ps(vector, _MM_SHUFFLE(3,3,3,3));
	x = _mm_mul_ps(x, matrix[0]);
	y = _mm_mul_ps(y, matrix[1]);
	z = _mm_mul_ps(z, matrix[2]);
	w = _mm_mul_ps(w, matrix[3]);

	return _mm_add_ps(_mm_add_ps(x, y), _mm_add_ps(z, w));
}


 sf_m128 CSIMD_SSE3::mat4x4_mul_sse(const sf_m128 *matrix, sf_m128 vector)
{
	sf_m128 x = _mm_mul_ps(matrix[0], vector);
	sf_m128 y = _mm_mul_ps(matrix[1], vector);
	sf_m128 z = _mm_mul_ps(matrix[2], vector);
	sf_m128 w = _mm_mul_ps(matrix[3], vector);
	sf_m128 tmp1 = _mm_hadd_ps(x, y); // = [y2+y3, y0+y1, x2+x3, x0+x1]
	sf_m128 tmp2 = _mm_hadd_ps(z, w); // = [w2+w3, w0+w1, z2+z3, z0+z1]

	return _mm_hadd_ps(tmp1, tmp2); // = [w0+w1+w2+w3, z0+z1+z2+z3, y0+y1+y2+y3, x0+x1+x2+x3]
}

sf_m128 CSIMD_SSE41::mat4x4_mul_sse(const sf_m128 *matrix, sf_m128 vector)
{
	sf_m128 x = _mm_dp_ps(matrix[0], vector, 0xF0 | 0x0F); // Choose to multiply x, y, z and w (0xF0 = 1111 0000), and store the output to all indices (0x0F == 0000 1111).
	sf_m128 y = _mm_dp_ps(matrix[1], vector, 0xF0 | 0x0F);
	sf_m128 z = _mm_dp_ps(matrix[2], vector, 0xF0 | 0x0F);
	sf_m128 w = _mm_dp_ps(matrix[3], vector, 0xF0 | 0x0F);

	sf_m128 xy = _mm_movelh_ps(x, y); // xy = [ _, y, _, x]
	sf_m128 zw = _mm_movelh_ps(z, w); // zw = [ _, w, _, z]

	return _mm_shuffle_ps(xy, zw, _MM_SHUFFLE(2, 0, 2, 0)); // ret = [w, z, y, x]
}

bool  VPCALL CSIMD_SSE::intersectLineAABB(const CAABBox &box, const CVec4D &rayPos, const CVec4D &rayDir, float tNear, float tFar)
{
	SMF_ASSERT(rayDir.isNormalized4());
	SMF_ASSERT(tNear <= tFar && "CAABBox::intersectLineAABB: User gave a degenerate line as input for the intersection test!");
	/* For reference, this is the C++ form of the vectorized SSE code below.

	CVec4D recipDir = rayDir.RecipFast4();
	CVec4D t1 = (aabbMinPoint - rayPos).mul(recipDir);
	CVec4D t2 = (aabbMaxPoint - rayPos).mul(recipDir);
	CVec4D near = t1.min(t2);
	CVec4D far = t1.MAX(t2);
	CVec4D rayDirAbs = rayDir.abs();

	if (rayDirAbs.x > 1e-4f) // ray is parallel to plane in question
	{
		tNear = MAX(near.x, tNear); // tNear tracks distance to intersect (enter) the CAABBox.
		tFar = min(far.x, tFar); // tFar tracks the distance to exit the CAABBox.
	}
	else if (rayPos.x < aabbMinPoint.x || rayPos.x > aabbMaxPoint.x) // early-out if the ray can't possibly enter the box.
		return false;

	if (rayDirAbs.y > 1e-4f) // ray is parallel to plane in question
	{
		tNear = MAX(near.y, tNear); // tNear tracks distance to intersect (enter) the CAABBox.
		tFar = min(far.y, tFar); // tFar tracks the distance to exit the CAABBox.
	}
	else if (rayPos.y < aabbMinPoint.y || rayPos.y > aabbMaxPoint.y) // early-out if the ray can't possibly enter the box.
		return false;

	if (rayDirAbs.z > 1e-4f) // ray is parallel to plane in question
	{
		tNear = MAX(near.z, tNear); // tNear tracks distance to intersect (enter) the CAABBox.
		tFar = min(far.z, tFar); // tFar tracks the distance to exit the CAABBox.
	}
	else if (rayPos.z < aabbMinPoint.z || rayPos.z > aabbMaxPoint.z) // early-out if the ray can't possibly enter the box.
		return false;

	return tNear < tFar;
	*/

	sf_m128 recipDir = _mm_rcp_ps(rayDir.v);
	// Note: The above performs an approximate reciprocal (11 bits of precision).
	// For a full precision reciprocal, perform a div:
//	sf_m128 recipDir = _mm_div_ps(_mm_set1_ps(1.f), rayDir.v);

	sf_m128 t1 = _mm_mul_ps(_mm_sub_ps(box.minPoint_SSE(), rayPos.v), recipDir);
	sf_m128 t2 = _mm_mul_ps(_mm_sub_ps(box.maxPoint_SSE(), rayPos.v), recipDir);

	sf_m128 nearD = _mm_min_ps(t1, t2); // [0 n3 n2 n1]
	sf_m128 farD = _mm_max_ps(t1, t2);  // [0 f3 f2 f1]

	// Check if the ray direction is parallel to any of the cardinal axes, and if so,
	// mask those [near, far] ranges away from the hit test computations.
	sf_m128 rayDirAbs = abs_ps(rayDir.v);

	const sf_m128 epsilon = _mm_set1_ps(1e-4f);
	// zeroDirections[i] will be nonzero for each axis i the ray is parallel to.
	sf_m128 zeroDirections = _mm_cmple_ps(rayDirAbs, epsilon);

	const sf_m128 floatInf = _mm_set1_ps(CMath::INFINITY_FLOAT);
	const sf_m128 floatNegInf = _mm_set1_ps(-CMath::INFINITY_FLOAT);

	// If the ray is parallel to one of the axes, replace the slab range for that axis
	// with [-inf, inf] range instead. (which is a no-op in the comparisons below)
	nearD = cmov_ps(nearD, floatNegInf, zeroDirections);
	farD = cmov_ps(farD , floatInf, zeroDirections);

	// Next, we need to compute horizontally max(nearD[0], nearD[1], nearD[2]) and min(farD[0], farD[1], farD[2])
	// to see if there is an overlap in the hit ranges.
	sf_m128 v1 = _mm_shuffle_ps(nearD, farD, _MM_SHUFFLE(0, 0, 0, 0)); // [f1 f1 n1 n1]
	sf_m128 v2 = _mm_shuffle_ps(nearD, farD, _MM_SHUFFLE(1, 1, 1, 1)); // [f2 f2 n2 n2]
	sf_m128 v3 = _mm_shuffle_ps(nearD, farD, _MM_SHUFFLE(2, 2, 2, 2)); // [f3 f3 n3 n3]
	nearD = _mm_max_ps(v1, _mm_max_ps(v2, v3));
	farD = _mm_min_ps(v1, _mm_min_ps(v2, v3));
	farD = _mm_shuffle_ps(farD, farD, _MM_SHUFFLE(3, 3, 3, 3)); // Unpack the result from high offset in the register.
	nearD = _mm_max_ps(nearD, _mm_set_ss(tNear));
	farD = _mm_min_ps(farD, _mm_set_ss(tFar));

	// Finally, test if the ranges overlap.
	sf_m128 rangeIntersects = _mm_cmple_ss(nearD, farD);

	// To store out out the interval of intersection, uncomment the following:
	// These are disabled, since without these, the whole function runs without a single memory store,
	// which has been profiled to be very fast! Uncommenting these causes an order-of-magnitude slowdown.
	// For now, using the SSE version only where the tNear and tFar ranges are not interesting.
//	_mm_store_ss(&tNear, nearD);
//	_mm_store_ss(&tFar, farD);

	// To avoid false positives, need to have an additional rejection test for each cardinal axis the ray direction
	// is parallel to.
	sf_m128 out2 = _mm_cmplt_ps(rayPos.v, box.minPoint_SSE());
	sf_m128 out3 = _mm_cmpgt_ps(rayPos.v, box.maxPoint_SSE());
	out2 = _mm_or_ps(out2, out3);
	zeroDirections = _mm_and_ps(zeroDirections, out2);

	sf_m128 yOut = _mm_shuffle_ps(zeroDirections, zeroDirections, _MM_SHUFFLE(1,1,1,1));
	sf_m128 zOut = _mm_shuffle_ps(zeroDirections, zeroDirections, _MM_SHUFFLE(2,2,2,2));

	zeroDirections = _mm_or_ps(_mm_or_ps(zeroDirections, yOut), zOut);
	// intersection occurs if the slab ranges had positive overlap and if the test was not rejected by the ray being
	// parallel to some cardinal axis.
	sf_m128 intersects = _mm_andnot_ps(zeroDirections, rangeIntersects);
	sf_m128 epsilonMasked = _mm_and_ps(epsilon, intersects);
	return _mm_comieq_ss(epsilon, epsilonMasked) != 0;
}

#define MATH_GEN_SSE2
#include "TriangleMesh_IntersectRay_SSE.inl"

#define MATH_GEN_SSE2
#define MATH_GEN_TRIANGLEINDEX
#include "TriangleMesh_IntersectRay_SSE.inl"

#define MATH_GEN_SSE2
#define MATH_GEN_TRIANGLEINDEX
#define MATH_GEN_UV
#include "TriangleMesh_IntersectRay_SSE.inl"

#undef MATH_GEN_SSE2
#undef MATH_GEN_TRIANGLEINDEX
#undef MATH_GEN_UV


#define MATH_GEN_SSE41
#include "TriangleMesh_IntersectRay_SSE.inl"

#define MATH_GEN_SSE41
#define MATH_GEN_TRIANGLEINDEX
#include "TriangleMesh_IntersectRay_SSE.inl"

#define MATH_GEN_SSE41
#define MATH_GEN_TRIANGLEINDEX
#define MATH_GEN_UV
#include "TriangleMesh_IntersectRay_SSE.inl"

#undef MATH_GEN_SSE2
#undef MATH_GEN_TRIANGLEINDEX
#undef MATH_GEN_UV


#define MATH_GEN_AVX
#include "TriangleMesh_IntersectRay_AVX.inl"

#define MATH_GEN_AVX
#define MATH_GEN_TRIANGLEINDEX
#include "TriangleMesh_IntersectRay_AVX.inl"

#define MATH_GEN_AVX
#define MATH_GEN_TRIANGLEINDEX
#define MATH_GEN_UV
#include "TriangleMesh_IntersectRay_AVX.inl"


} //end MATH

} //end SMF
