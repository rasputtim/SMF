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
#include <xmmintrin.h>
__attribute__((  aligned(16)  ))  float _SIMDx86_float_SSE_NO_W_MASK[4]={0xffffffff,0xffffffff,0xffffffff,0x0,};

#define HIPREC
//-----------------------------------------------------------------------------
// SSE implementations of optimized routines:
//-----------------------------------------------------------------------------
float _SSE_Sqrt(float x)
{
	float	root = 0.f;

	__asm__ __volatile__(
		"movss %1,%%xmm2\n"
		"sqrtss %%xmm2,%%xmm1\n"
		//"movss %%xmm1,%0"
       	"movss %%xmm2, -4(%%esp)\n"
		"flds -4(%%esp)\n": "=t" (root)
		: "m" (x)
	);

	return root;
}

// Single iteration NewtonRaphson reciprocal square root:
// 0.5 * rsqrtps * (3 - x * rsqrtps(x) * rsqrtps(x))
// Very low error, and fine to use in place of 1.f / sqrtf(x).

// Intel / Kipps SSE RSqrt.  Significantly faster than above.
float _SSE_RSqrtAccurate(float a)
{
	__attribute__((  aligned(16)  )) float half = 0.5f;
	__attribute__((  aligned(16)  )) float three = 3.f;
    __attribute__((  aligned(16)  )) float x;

#if defined (MSCVER)
	__asm
	{
		movss   xmm3, a;
		movss   xmm1, half;
		movss   xmm2, three;
		rsqrtss xmm0, xmm3;

		mulss   xmm3, xmm0;
		mulss   xmm1, xmm0;
		mulss   xmm3, xmm0;
		subss   xmm2, xmm3;
		mulss   xmm1, xmm2;

		movss   x,    xmm1;
	}
#endif // _WIN32
#if !defined(MASM_INTEL)
	__asm__ __volatile__(
		"movss   %1, %%xmm3 \n"
        "movss   %2, %%xmm1 \n"
        "movss   %3, %%xmm2 \n"
        "rsqrtss %%xmm3, %%xmm0 \n"
        "mulss   %%xmm0, %%xmm3 \n"
        "mulss   %%xmm0, %%xmm1 \n"
        "mulss   %%xmm0, %%xmm3 \n"
        "subss   %%xmm3, %%xmm2 \n"
        "mulss   %%xmm2, %%xmm1 \n"
        "movss   %%xmm1, %0 \n"
		: "=m" (x)
		: "m" (a), "m" (half), "m" (three)
);
#else  //masm=intel
	#error "Not Implemented"
#endif

	return x;
}


// Simple SSE rsqrt.  Usually accurate to around 6 (relative) decimal places
// or so, so ok for closed transforms.  (ie, computing lighting normals)
float _SSE_RSqrtFast(float x)
{

	__attribute__((  aligned(16)  )) float rroot;
#if defined(MSCVER)
	_asm
	{
		rsqrtss	xmm0, x
		movss	rroot, xmm0
	}
#endif // _WIN32
#if !defined(MASM_INTEL)
	 asm(
		"rsqrtss %1, %%xmm0 \n"
		//"rsqrtss %%xmm1, %%xmm1\n"
		//"movss %%xmm0, %0 \n\t"
		"movss %%xmm0, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		: "=t" (x)
		: "m" (rroot)
		:
	);
#else  //masm_intel

#endif

	return rroot;
}
inline float
sqrtf ( float x )
{
#ifdef PURE_VECTOR
    __m128  xx = _mm_load_ss( & x );
    __m128  xr = _mm_rsqrt_ss( xx );
    __m128  xt;

    xt = _mm_mul_ss( xr, xr );
    xt = _mm_mul_ss( xt, xx );
    xt = _mm_sub_ss( f3, xt );
    xt = _mm_mul_ss( xt, f05 );
    xr = _mm_mul_ss( xr, xt );

    _mm_store_ss( & x, xr );

    return x;
#else
    float   r;

    _mm_store_ss( & r, _mm_sqrt_ss( _mm_load_ss( & x ) ) );

    r *= ((3.0f - r * r * x) * 0.5f);

    return r;
#endif
}
float SqRt_SSE( const float* x )
    {
        __attribute__((  aligned(16)  ))  __m128 r0;
        __attribute__((  aligned(16)  ))  float V[4]={14.0f,0.0f,0.0f,0.0f};
        __attribute__((  aligned(16)  ))  __m128 v_i;
        //v_i = _mm_loadu_ps(x);
        __attribute__((  aligned(16)  )) float teste=14.f;
        r0 = _mm_load_ss( &teste );
        r0 = _mm_sqrt_ss( r0 );

        float y;

        _mm_store_ss( &y, r0 );

        return y;
    }
float  Vector3D_Length(const float* pVec)
{




		__attribute__((  aligned(16)  )) float dummy;
		/* SSE/SSE2 Implementation */
		asm(
		"movaps (%1), %%xmm0\n"
		"andps %2, %%xmm0\n"    /* xmm0 = x | y | z | 0.0f */
		"mulps %%xmm0, %%xmm0\n"

		/* Shift data around in the registers (I loves me some haddps right now!!) */
		"movhlps %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   |  z's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | 0   | z's + x's */
		"shufps $0x55, %%xmm0, %%xmm0\n"/* xmm0 = y's | y's | y's | y's */
		"addss %%xmm0, %%xmm1\n"		/* xmm1 = ?   | ?   | ?   | x's+y's+z's */

		#ifdef HIPREC
		/* Full square root */
		"sqrtss %%xmm1, %%xmm1\n"
		#else
		/* rcp( rsqrt(value) ) (This may be very inaccurate) */
		//"rsqrtss %%xmm1, %%xmm1\n"
		//"rcpss %%xmm1, %%xmm1\n"
		#endif

		"movss %%xmm1, -4(%%esp)\n"
		"flds -4(%%esp)\n"
		: "=t" (dummy)
		: "r" (pVec),"m"(_SIMDx86_float_SSE_NO_W_MASK)
		);
        //float test = __builtin_sqrt (dummy);
		#ifdef HIPREC
//		 __builtin_ia32_sqrtss(raiz);
		//return _SSE_Sqrt(dummy);
        return dummy;
        #else
        return _SSE_RSqrtAccurate(dummy);
        //return dummy;
        #endif
}
//This is for the normal PC
//The programm entry point
int main(int argc,char *argv[])
{
__attribute__((  aligned(16)  )) float var1[3];
__attribute__((  aligned(16)  )) float var2[4];
__attribute__((  aligned(16)  )) float var3[4];
__attribute__((  aligned(16)  )) float result;
var1[0]=1.0f;
var1[1]=2.0f;
var1[2]=3.0f;

float * varPTR = (float *) &var1;

var3[0]=1.0f;
var3[1]=1.0f;
var3[2]=1.0f;
var3[3]=1.0f;
int count=4;

printf("Before Lenght -> %f %f %f\n", var1[0], var1[1], var1[2]);
result = Vector3D_Length(var1);
printf("After Length -> %f %f %f \n", var1[0], var1[1], var1[2]);
printf("Lenght is -> %f \n", result);

return 0;
}



