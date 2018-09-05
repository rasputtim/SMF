//    Open Mugen is a redevelopment of Elecbyte's M.U.G.E.N wich will be 100% compatible to it
//    Copyright (C) 2004  Sahin Vardar
//
//    If you know bugs or have a wish on Open Muegn or (money/girls/a car) for me ;-)
//    Feel free and email me: sahin_v@hotmail.com  ICQ:317502935
//    Web: http://openmugen.sourceforge.net/
//    --------------------------------------------------------------------------
//
//    This program is free software; you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation; either version 2 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.

//    You should have received a copy of the GNU General Public License
//    along with this program; if not, write to the Free Software
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.



#include <io.h>
#include <string>
#include <vector>
#include <queue>
#include <list>
#include <iostream>

void Copy16( float *dst, const float *src, const int count ) {
//	Debug::debug(Debug::math, __FUNCTION__) << "size of float: " << sizeof(float) << endl;

	__asm__ (
		"mov		ecx, %0\n"
		"mov		edx, %1\n"
		"mov		eax, %2\n"
		"add		eax, 3\n"
		"shr		eax, 2\n"
		"jz		doneCopy16\n"
		"shl		eax, 4\n"
		"add		ecx, eax\n"
		"add		edx, eax\n"
		"neg		eax\n"
	"loopCopy16:\n"
		"movaps	xmm0, [ecx+eax]\n"
		"movaps	[edx+eax], xmm0\n"
		"add		eax, 16\n"
		"jl		loopCopy16\n"
	"doneCopy16:\n"
	:
	    : "m" (src), "m" (dst),"m"(count)
	        :);
}

#define ALIGNED  __attribute__((  aligned(16)  ))
#define ALIGN4_INIT1( X, INIT )				ALIGNED static X[4] = { INIT, INIT, INIT, INIT }
#define ALIGN4_INIT4( X, I0, I1, I2, I3 )	ALIGNED static X[4] = { I0, I1, I2, I3 }
#define ALIGN2_INIT1( X, INIT )				ALIGNED static X[2] = { INIT, INIT }
#define ALIGN2_INIT2( X, I0, I1 )			ALIGNED static X[2] = { I0, I1 }

const float	PI				= 3.14159265358979323846f;
const float	TWO_PI			= 2.0f * PI;
const float	HALF_PI			= 0.5f * PI;
const float	ONEFOURTH_PI	= 0.25f * PI;
const float E				= 2.71828182845904523536f;
const float SQRT_TWO		= 1.41421356237309504880f;
const float SQRT_THREE		= 1.73205080756887729352f;
const float	SQRT_1OVER2		= 0.70710678118654752440f;
const float	SQRT_1OVER3		= 0.57735026918962576450f;
const float	M_DEG2RAD		= PI / 180.0f;
const float	M_RAD2DEG		= 180.0f / PI;
const float	M_SEC2MS		= 1000.0f;
const float	M_MS2SEC		= 0.001f;
const float	INFINITY		= 1e30f;
const float FLT_EPSILON		= 1.192092896e-07f;


ALIGN4_INIT4( unsigned long SIMD_SP_singleSignBitMask, (unsigned long) ( 1 << 31 ), 0, 0, 0 );
ALIGN4_INIT1( unsigned long SIMD_SP_signBitMask, (unsigned long) ( 1 << 31 ) );
ALIGN4_INIT1( unsigned long SIMD_SP_absMask, (unsigned long) ~( 1 << 31 ) );
ALIGN4_INIT1( unsigned long SIMD_SP_infinityMask, (unsigned long) ~( 1 << 23 ) );
ALIGN4_INIT1( unsigned long SIMD_SP_not, 0xFFFFFFFF );

ALIGN4_INIT1( float SIMD_SP_zero, 0.0f );
ALIGN4_INIT1( float SIMD_SP_half, 0.5f );
ALIGN4_INIT1( float SIMD_SP_one, 1.0f );
ALIGN4_INIT1( float SIMD_SP_two, 2.0f );
ALIGN4_INIT1( float SIMD_SP_three, 3.0f );
ALIGN4_INIT1( float SIMD_SP_four, 4.0f );
ALIGN4_INIT1( float SIMD_SP_maxShort, (1<<15) );
ALIGN4_INIT1( float SIMD_SP_tiny, 1e-10f );
ALIGN4_INIT1( float SIMD_SP_PI, PI );
ALIGN4_INIT1( float SIMD_SP_halfPI, HALF_PI );
ALIGN4_INIT1( float SIMD_SP_twoPI, TWO_PI );
ALIGN4_INIT1( float SIMD_SP_oneOverTwoPI, 1.0f / TWO_PI );
ALIGN4_INIT1( float SIMD_SP_infinity, INFINITY );
ALIGN4_INIT4( float SIMD_SP_lastOne, 0.0f, 0.0f, 0.0f, 1.0f );

ALIGN4_INIT1( float SIMD_SP_rsqrt_c0,  3.0f );
ALIGN4_INIT1( float SIMD_SP_rsqrt_c1, -0.5f );
ALIGN4_INIT1( float SIMD_SP_mat2quat_rsqrt_c1, -0.5f*0.5f );

ALIGN4_INIT1( float SIMD_SP_sin_c0, -2.39e-08f );
ALIGN4_INIT1( float SIMD_SP_sin_c1,  2.7526e-06f );
ALIGN4_INIT1( float SIMD_SP_sin_c2, -1.98409e-04f );
ALIGN4_INIT1( float SIMD_SP_sin_c3,  8.3333315e-03f );
ALIGN4_INIT1( float SIMD_SP_sin_c4, -1.666666664e-01f );

ALIGN4_INIT1( float SIMD_SP_cos_c0, -2.605e-07f );
ALIGN4_INIT1( float SIMD_SP_cos_c1,  2.47609e-05f );
ALIGN4_INIT1( float SIMD_SP_cos_c2, -1.3888397e-03f );
ALIGN4_INIT1( float SIMD_SP_cos_c3,  4.16666418e-02f );
ALIGN4_INIT1( float SIMD_SP_cos_c4, -4.999999963e-01f );

ALIGN4_INIT1( float SIMD_SP_atan_c0,  0.0028662257f );
ALIGN4_INIT1( float SIMD_SP_atan_c1, -0.0161657367f );
ALIGN4_INIT1( float SIMD_SP_atan_c2,  0.0429096138f );
ALIGN4_INIT1( float SIMD_SP_atan_c3, -0.0752896400f );
ALIGN4_INIT1( float SIMD_SP_atan_c4,  0.1065626393f );
ALIGN4_INIT1( float SIMD_SP_atan_c5, -0.1420889944f );
ALIGN4_INIT1( float SIMD_SP_atan_c6,  0.1999355085f );
ALIGN4_INIT1( float SIMD_SP_atan_c7, -0.3333314528f );

ALIGN4_INIT1(double _SIMDx86_double_one, 1.0);
ALIGN4_INIT1(float  _SIMDx86_float_one, 1.0f);
ALIGN4_INIT4(float  _SIMDx86_float_one_w, 0.0f, 0.0f, 0.0f, 1.0f );
ALIGN4_INIT4(float  _SIMDx86_float_one_y, 0.0f, 1.0f, 0.0f, 0.0f );

ALIGN2_INIT2(unsigned int  _SIMDx86_float_POSNEG, 0x00000000, 0x80000000 );	/* XOR Changes the signs to + - */
ALIGN2_INIT2(unsigned int  _SIMDx86_float_NEGPOS, 0x80000000, 0x00000000 );	/* XOR Changes the signs to - + */
ALIGN2_INIT2(unsigned int  _SIMDx86_float_NEGNEG, 0x80000000, 0x80000000 );	/* XOR Changes the signs to - - */
ALIGN2_INIT2(unsigned int  _SIMDx86_float_3DNOW_NO_W_MASK, 0xFFFFFFFF, 0x00000000 );
ALIGN4_INIT4(unsigned int  _SIMDx86_float_ABS, 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF,0x7FFFFFFF );


ALIGNED  float _SIMDx86_float_SSE_NO_W_MASK[4]={0xffffffff,0xffffffff,0xffffffff,0x0};

float  Vector4D_Dot(const float* pSrc1, const float* pSrc2)
{

	/* SSE/SSE2 Implementation */
	float dummy;
	float *min=&dummy;
		asm(
			"mov eax, %0\n"                       // Load pointers into CPU regs and multiply
			"mov ebx, %1\n"
			"movups xmm0, [eax]\n"
			"movups xmm1, [ebx]\n"
			"mulps xmm1, xmm0\n"
			"movaps xmm2,xmm1\n"
			"shufps xmm2, xmm2,  0x1B\n" /*0x1B= 00011011 / xmm2 = x | y | z | w */

			"addps  xmm2, xmm1\n"			/* xmm2 = w+x | y+z | y+z | w+x */
			"movss  xmm3, xmm2\n"			/* xmm3 = ??? | ??? | ??? | w+x */
			"shufps xmm1, xmm1, 0x01\n"  /* 0x01= 00000001 / xmm3 = ??? | ??? | ??? | w+x */
			"addss  xmm2, xmm3\n"			/* xmm2 = ??? | ??? | ??? | dot4 */


			//coloca 0-31 de xmm2 em dummy
			"mov			esi, %2\n"
			"movss		[esi], xmm2\n"

	:
	: "r" (pSrc1), "r" (pSrc2),"m" (min)
	);

	return dummy;



}

void Vector3D_AlignedSumOf(float* pOut, const float* pIn1, const float* pIn2)
{
		/* SSE/SSE2/SSE3 Implementation */
		asm(
                "mov eax, %0\n"                       // Load pointers into CPU regs
                "mov ebx, %1\n"
				"mov ecx, %2\n"
                "movaps xmm0, [eax]\n"                 // Move unaligned vectors to SSE regs
                "movaps xmm1, [ebx]\n"
				"andps xmm0, %3\n" //zera o elemento w
                "andps xmm1, %3\n"  //re"move w component
				"addps xmm0, xmm1\n"                   // Add vector elements
                //"movlpd	[ecx+ 0], xmm0\n"  //retira os tres bytes
				//"shufps xmm0,xmm0, 0x02\n" //shuffle colocar nos bits 0-31 o 3o byte
				//"movss	[ecx+ 8], xmm0\n"
                /* Store*/
                "movlpd [%2] , xmm0\n"
                "shufps  xmm0, xmm0, 0x02\n"  //shuffle colocar nos bits 0-31 o 3o byte
                "movss	[%2+ 8], xmm0\n"
	:
	: "r" (pIn1), "r" (pIn2), "r" (pOut), "m" (_SIMDx86_float_SSE_NO_W_MASK)
	);


}

float SSE_Sin( float a ) {
#if 1

	float t;

	asm (
		"movss		xmm1, %13\n"
		"movss		xmm2, xmm1\n"
		"movss		xmm3, xmm1\n"
		"mulss		xmm2, %0\n"
		"cvttss2si	ecx, xmm2\n"
		"cmpltss	xmm3, %1\n"
		"andps		xmm3, %2\n"
		"cvtsi2ss	xmm2, ecx\n"
		"subss		xmm2, xmm3\n"
		"mulss		xmm2, %3\n"
		"subss		xmm1, xmm2\n"

		"movss		xmm0, %4\n"			// xmm0 = PI
		"subss		xmm0, xmm1\n"					// xmm0 = PI - a
		"movss		xmm1, xmm0\n"					// xmm1 = PI - a
		"andps		xmm1, %5\n"	// xmm1 = signbit( PI - a )
		"movss		xmm2, xmm0\n"					// xmm2 = PI - a
		"xorps		xmm2, xmm1\n"					// xmm2 = fabs( PI - a )
		"cmpnltss	xmm2, %6\n"		// xmm2 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? 0xFFFFFFFF : 0x00000000
		"movss		xmm3, %4\n"			// xmm3 = PI
		"xorps		xmm3, xmm1\n"					// xmm3 = PI ^ signbit( PI - a )
		"andps		xmm3, xmm2\n"					// xmm3 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? ( PI ^ signbit( PI - a ) ) : 0.0f
		"andps		xmm2, %5\n"	// xmm2 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? SIMD_SP_signBitMask : 0.0f
		"xorps		xmm0, xmm2\n"
		"addps		xmm0, xmm3\n"

		"movss		xmm1, xmm0\n"
		"mulss		xmm1, xmm1\n"
		"movss		xmm2, %7\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %8\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %9\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %10\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %11\n"
		"mulss		xmm2, xmm1\n"
		"addss		xmm2, %12\n"
		"mulss		xmm2, xmm0\n"
		/* Store */
		"movss		%14, xmm2\n"
	    :
        : "m" (SIMD_SP_oneOverTwoPI), "m"(SIMD_SP_zero),"m"(SIMD_SP_one),"m"(SIMD_SP_twoPI),
          "m"(SIMD_SP_PI),"m"(SIMD_SP_signBitMask), "m"(SIMD_SP_halfPI),
          "m"(SIMD_SP_sin_c0),"m"(SIMD_SP_sin_c1),"m"(SIMD_SP_sin_c2),
          "m"(SIMD_SP_sin_c3),"m"(SIMD_SP_sin_c4),"m"(SIMD_SP_one),"m"(a),"m"(t)
        );
	return t;

#else

	float s, t;

	if ( ( a < 0.0f ) || ( a >= CMath::TWO_PI ) ) {
		a -= floorf( a / CMath::TWO_PI ) * CMath::TWO_PI;
	}

	a = CMath::PI - a;
	if ( fabs( a ) >= CMath::HALF_PI ) {
		a = ( ( a < 0.0f ) ? -CMath::PI : CMath::PI ) - a;
	}

	s = a * a;
	t = -2.39e-08f;
	t *= s;
	t += 2.7526e-06f;
	t *= s;
	t += -1.98409e-04f;
	t *= s;
	t += 8.3333315e-03f;
	t *= s;
	t += -1.666666664e-01f;
	t *= s;
	t += 1.0f;
	t *= a;

	return t;
#endif

}

void SSE_Sin4( float entrada[4], float saida[4] ) {

    asm (
		"mov		edi, %13\n"
		"mov		esi, %14\n"
		"movaps		xmm1, [edi]\n"
		"movaps		xmm2, xmm1\n"
		"mulps		xmm2, %0\n"
		"movhlps	xmm3, xmm2\n"
		"cvttss2si	ecx, xmm2\n"
		"cvtsi2ss	xmm2, ecx\n"
		"cvttss2si	edx, xmm3\n"
		"cvtsi2ss	xmm3, edx\n"
		"shufps		xmm2, xmm2, 0x01\n"
		"shufps		xmm3, xmm3, 0x01\n"
		"cvttss2si	ecx, xmm2\n"
		"cvtsi2ss	xmm2, ecx\n"
		"cvttss2si	edx, xmm3\n"
		"cvtsi2ss	xmm3, edx\n"
		"shufps		xmm2, xmm3,0x11 \n"   /*R_SHUFFLEPS( 1, 0, 1, 0 )*/
		"movaps		xmm3, xmm1\n"
		"cmpltps	xmm3, %1\n"
		"andps		xmm3, %2\n"
		"subps		xmm2, xmm3\n"
		"mulps		xmm2, %3\n"
		"subps		xmm1, xmm2\n"

		"movaps		xmm0, %4\n"			// xmm0 = PI
		"subps		xmm0, xmm1\n"					// xmm0 = PI - a
		"movaps		xmm1, xmm0\n"					// xmm1 = PI - a
		"andps		xmm1, %5\n"	// xmm1 = signbit( PI - a )
		"movaps		xmm2, xmm0\n"					// xmm2 = PI - a
		"xorps		xmm2, xmm1\n"					// xmm2 = fabs( PI - a )
		"cmpnltps	xmm2, %6\n"	// xmm2 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? 0xFFFFFFFF : 0x00000000
		"movaps		xmm3, %4\n"			// xmm3 = PI
		"xorps		xmm3, xmm1\n"					// xmm3 = PI ^ signbit( PI - a )
		"andps		xmm3, xmm2\n"					// xmm3 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? ( PI ^ signbit( PI - a ) ) : 0.0f
		"andps		xmm2, %5\n"	// xmm2 = ( fabs( PI - a ) >= CMath::HALF_PI ) ? SIMD_SP_signBitMask : 0.0f
		"xorps		xmm0, xmm2\n"
		"addps		xmm0, xmm3\n"

		"movaps		xmm1, xmm0\n"
		"mulps		xmm1, xmm1\n"
		"movaps		xmm2, %7\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %8\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %9\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %10\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %11\n"
		"mulps		xmm2, xmm1\n"
		"addps		xmm2, %12\n"
		"mulps		xmm2, xmm0\n"
		"movaps		[esi], xmm2\n"
        :
	    : "m" (SIMD_SP_oneOverTwoPI), "m"(SIMD_SP_zero),"m"(SIMD_SP_one),"m"(SIMD_SP_twoPI),
          "m"(SIMD_SP_PI),"m"(SIMD_SP_signBitMask), "m"(SIMD_SP_halfPI),
          "m"(SIMD_SP_sin_c0),"m"(SIMD_SP_sin_c1),"m"(SIMD_SP_sin_c2),
          "m"(SIMD_SP_sin_c3),"m"(SIMD_SP_sin_c4),"m"(SIMD_SP_one),"r"(entrada),"r"(saida)
        );

}
//This is for the normal PC
//The programm entry point
int main(int argc,char *argv[])
{
__attribute__((  aligned(16)  )) float var1[4];
__attribute__((  aligned(16)  )) float var2[4];
__attribute__((  aligned(16)  )) float var3[4];
var1[0]=2.0f;
var1[1]=3.0f;
var1[2]=4.0f;
var1[3]=5.0f;
var3[0]=1.0f;
var3[1]=1.0f;
var3[2]=1.0f;
var3[3]=1.0f;
int count=4;
float sin;
printf("Before Copy-> %f %f %f %f\n", var2[0], var2[1], var2[2], var2[3]);
Copy16(var2,var1,count);

printf("After Copy -> %f %f %f %f\n", var2[0], var2[1], var2[2], var2[3]);
Vector3D_AlignedSumOf(var3,var1,var2);
printf("After Sum -> %f %f %f %f\n", var3[0], var3[1], var3[2], var3[3]);
float x=Vector4D_Dot(var3,var1);
printf("After Dot -> %f %f %f %f\n", var3[0], var3[1], var3[2], var3[3]);
printf("Dot is -> %f \n", x);
var3[0]=1.0f;
var3[1]=1.0f;
var3[2]=1.0f;
var3[3]=1.0f;

SSE_Sin4( var3,var2 );
printf("Sin -> %f %f %f %f\n", var2[0], var2[1], var2[2], var2[3]);
sin=SSE_Sin(1.0f);
printf("Sin -> %f\n", sin);


    return 0;

}



