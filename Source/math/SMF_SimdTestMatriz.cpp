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

#include "math/SMF_Simd.h"
//#pragma hdrstop

#include "math/SMF_SimdGeneric.h"
#include "math/SMF_SimdMMX.h"
#include "math/SMF_Simd3DNow.h"
#include "math/SMF_SimdSSE.h"
#include "math/SMF_SimdSSE2.h"
#include "math/SMF_SimdSSE3.h"
//#include "math/SMF_SimdSSE41.h"
//#include "math/SMF_SimdSSE42.h"
//#include "math/SMF_SimdAltiVec.h"
#include "util/SMF_Random.h"
#include "math/SMF_Math.h"
#include "math/SMF_Vector.h"
#include "math/SMF_Matriz.h"
#include "geometry/SMF_Plane.h"
#include "math/SMF_EulerAngles.h"
#include "util/SMF_StringUtils.h"
#include "geometry/SMF_DrawVert.h"
#include "math/SMF_JointTransform.h"
#include "sys/SMF_System.h"
#include "util/SMF_Debug.h"

namespace SMF {
namespace MATH{

//===============================================================
//
// Test code
//
//===============================================================

#define COUNT		1024		// data count
#define NUMTESTS	2048		// number of tests
#define	MATX_SIMD_EPSILON 1e-5f

#define RANDOM_SEED		1013904223L	//((int)idLib::sys->GetClockTicks())
extern CSIMDProcessor *p_simd;
extern CSIMDProcessor *p_generic;
extern CSIMDProcessor *	processor;			// pointer to SIMD processor
extern CSIMDProcessor *	generic;				// pointer to generic SIMD implementation
extern CSIMDProcessor *	SIMDProcessor;
extern unsigned baseClocks;

//unsigned saved_ebx = 0;
//unsigned start_ClockCount = 0;
//unsigned end_ClockCount = 0;
//double ticksPerNanosecond;

typedef struct {
	unsigned int hi;
	unsigned int lo;
} U64;


#if defined(_WIN32) && defined(_MSC_VER)

#define TIME_TYPE int

#pragma warning(disable : 4731)     // frame pointer register 'ebx' modified by inline assembly code




#elif MACOS_X

#include <stdlib.h>
#include <unistd.h>			// this is for sleep()
#include <sys/time.h>
#include <sys/resource.h>
#include <mach/mach_time.h>

double ticksPerNanosecond;

#define TIME_TYPE uint64_t

#ifdef __MWERKS__ //time_in_millisec is missing
/*

    .text
	.align 2
	.globl _GetTB
_GetTB:

loop:
	        mftbu   r4	;  load from TBU
	        mftb    r5	;  load from TBL
	        mftbu   r6	;  load from TBU
	        cmpw    r6, r4	;  see if old == new
	        bne     loop	;  if not, carry occured, therefore loop

	        stw     r4, 0(r3)
	        stw     r5, 4(r3)

done:
	        blr		;  return

*/


asm void GetTB(U64 *in)
{
	nofralloc			// suppress prolog
	machine 603			// allows the use of mftb & mftbu functions

loop:
	mftbu	r5			// grab the upper time base register (TBU)
	mftb	r4			// grab the lower time base register (TBL)
	mftbu	r6			// grab the upper time base register (TBU) again

	cmpw	r6,r5		// see if old TBU == new TBU
	bne-	loop		// loop if carry occurred (predict branch not taken)

	stw  	r4,4(r3)	// store TBL in the low 32 bits of the return value
	stw  	r5,0(r3)	// store TBU in the high 32 bits of the return value

	blr
}




double TBToDoubleNano( U64 startTime, U64 stopTime, double ticksPerNanosecond );

#if __MWERKS__
asm void GetTB( U64 * );
#else
void GetTB( U64 * );
#endif

double TBToDoubleNano( U64 startTime, U64 stopTime, double ticksPerNanosecond ) {
	#define K_2POWER32 4294967296.0
	#define TICKS_PER_NANOSECOND 0.025
	double nanoTime;
	U64 diffTime;

	// calc the difference in TB ticks
	diffTime.hi = stopTime.hi - startTime.hi;
	diffTime.lo = stopTime.lo - startTime.lo;

	// convert TB ticks into time
	nanoTime = (double)(diffTime.hi)*((double)K_2POWER32) + (double)(diffTime.lo);
	nanoTime = nanoTime/ticksPerNanosecond;
	return (nanoTime);
}

TIME_TYPE time_in_millisec() {
	#define K_2POWER32 4294967296.0
	#define TICKS_PER_NANOSECOND 0.025

	U64 the_time;
	double nanoTime, milliTime;

	GetTB( &the_time );

	// convert TB ticks into time
	nanoTime = (double)(the_time.hi)*((double)K_2POWER32) + (double)(the_time.lo);
	nanoTime = nanoTime/ticksPerNanosecond;

	// nanoseconds are 1 billionth of a second. I want milliseconds
	milliTime = nanoTime * 1000000.0;

	printf( "ticks per nanosec -- %lf\n", ticksPerNanosecond );
	printf( "nanoTime is %lf -- milliTime is %lf -- as int is %i\n", nanoTime, milliTime, (int)milliTime );

	return (int)milliTime;
}

#define StartRecordTimeLocal( start )			\
	start = time_in_millisec();

#define StopRecordTimeLocal( end )				\
	end = time_in_millisec();


#else
#define StartRecordTimeLocal( start )			\
	start = mach_absolute_time();

#define StopRecordTimeLocal( end )				\
	end = mach_absolute_time();
#endif


#else   //mingw version

#define TIME_TYPE unsigned int


#define GetBest( start, end, best )			\
	if ( !best || end - start < best ) {	\
		best = end - start;					\
	}

#endif

//===================CMat4D====================================

CVec3D outVec, inVec1, inVec2;
ALIGNTO16 CVec3D outVecAligned;
ALIGNTO16 CVec3D inVec1Aligned;
ALIGNTO16 CVec3D inVec2Aligned;

static void reset3D()
{
	outVec.toZero();
	outVecAligned.toZero();

	inVec1.x = 10.0f;
	inVec1.y = 20.0f;
	inVec1.z = 30.0f;

	inVec2.x = 1.0f;
	inVec2.y = 2.0f;
	inVec2.z = 3.0f;

	inVec1Aligned.x = 10.0f;
	inVec1Aligned.y = 20.0f;
	inVec1Aligned.z = 30.0f;

	inVec2Aligned.x = 1.0f;
	inVec2Aligned.y = 2.0f;
	inVec2Aligned.z = 3.0f;

}

CVec4D outVec4d, inVec4d1, inVec4d2;
ALIGNTO16  CVec4D outVec4dAligned;
ALIGNTO16  CVec4D inVec4d1Aligned;
ALIGNTO16  CVec4D inVec4d2Aligned;
static void reset4D()
{
	outVec4d.x = outVec4d.y = outVec4d.z= outVec4d.w = 0.0f;

	inVec4d1.x = 10.0f;
	inVec4d1.y = 20.0f;
	inVec4d1.z = 30.0f;
	inVec4d1.w = 40.0f;

	inVec4d2.x = 1.0f;
	inVec4d2.y = 2.0f;
	inVec4d2.z = 3.0f;
	inVec4d2.w = 4.0f;

	inVec4d1Aligned.x = 10.0f;
	inVec4d1Aligned.y = 20.0f;
	inVec4d1Aligned.z = 30.0f;
	inVec4d1Aligned.w = 40.0f;

	inVec4d2Aligned.x = 1.0f;
	inVec4d2Aligned.y = 2.0f;
	inVec4d2Aligned.z = 3.0f;
	inVec4d2Aligned.w = 4.0f;

}


CMat4D matout4d, matin4d1, matin4d2;
ALIGNTO16  CMat4D matout4dAligned;
ALIGNTO16  CMat4D matin4d1Aligned;
ALIGNTO16  CMat4D matin4d2Aligned;
static void matReset4D()
{
	matout4d.toZero();
	matout4dAligned.toZero();
	CVec4D vec1(-2,3,1,-1);
	CVec4D vec2(0,1,2,3);
	CVec4D vec3(1,-1,1,-2);
	CVec4D vec4(4,-3,5,1);

	CVec4D vec5(1,0,2,-1);
	CVec4D vec6(2,1,3,-2);
	CVec4D vec7(0,0,2,3);
	CVec4D vec8(1,-1,0,2);

	CMat4D mat1(vec1,vec2,vec3,vec4);
	CMat4D mat2(vec5,vec6,vec7,vec8);
	matin4d1 = mat1;
	matin4d2 = mat2;

	matin4d1Aligned = mat1;
	matin4d2Aligned = mat2;

}


static void TestMatriz4DSum()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CMat4D tst,dst;
	const char *result;
	matin4d1.toZero();
	matin4d2.toZero();
	matReset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D Sum()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->mat4D_Sum(&matin4d1, &matin4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = matin4d1;

		CProcClock::printClocks( va( "generic->Matriz4D_Sum "), 1, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->mat4D_Sum(&matin4d1, &matin4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = matin4d1;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->Matriz4D_Sum  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestMatriz4DSumOf()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CMat4D tst,dst;
	const char *result;
	matin4d1.toZero();
	matin4d2.toZero();
	matReset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D SumOf()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->mat4D_SumOf(&matout4d,&matin4d1, &matin4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = matin4d1;

		CProcClock::printClocks( va( "generic->Matriz4D_SumOf "), 1, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->mat4D_SumOf(&matout4d,&matin4d1, &matin4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = matin4d1;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->Matriz4D_SumOf  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestMatriz4DDiff()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CMat4D tst,dst;
	const char *result;
	matin4d1.toZero();
	matin4d2.toZero();
	matReset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D Diff()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->mat4D_Diff(&matin4d1, &matin4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = matin4d1;

		CProcClock::printClocks( va( "generic->Matriz4D_Diff "), 1, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->mat4D_Diff(&matin4d1, &matin4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = matin4d1;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->Matriz4D_Diff  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestMatriz4DDiffOf()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CMat4D tst,dst;
	const char *result;
	matin4d1.toZero();
	matin4d2.toZero();
	matReset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D DiffOf()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->mat4D_DiffOf(&matout4d,&matin4d1, &matin4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = matout4d;

		CProcClock::printClocks( va( "generic->Matriz4D_DiffOf "), 1, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->mat4D_DiffOf(&matout4d,&matin4d1, &matin4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = matout4d;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->Matriz4D_DiffOf  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestMatriz4DMultiply()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CMat4D tst,dst;
	const char *result;
	matin4d1.toZero();
	matin4d2.toZero();
	matReset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D multiply()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->mat4D_Multiply(&matin4d1, &matin4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = matin4d1;

		CProcClock::printClocks( va( "generic->Matriz4D_Multiply "), 1, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->mat4D_Multiply(&matin4d1, &matin4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = matin4d1;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->Matriz4D_Multiply  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestMatriz4DMultiplyOf()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CMat4D tst,dst;
	const char *result;
	matin4d1Aligned.toZero();
	matin4d2Aligned.toZero();
	matReset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D MultiplyOf()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->mat4D_MultiplyOf(&matout4d,&matin4d1, &matin4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = matout4d;

		CProcClock::printClocks( va( "generic->Matriz4D_MultiplyOf "), 1, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->mat4D_MultiplyOf(&matout4d,&matin4d1, &matin4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = matout4d;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->Matriz4D_MultiplyOf  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}


static void TestMatriz4DTranspose()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CMat4D tst,dst;
	const char *result;
	matin4d1.toZero();
	matin4d2.toZero();
	matReset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D transpose()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->mat4D_Transpose(&matin4d1);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = matin4d1;

		CProcClock::printClocks( va( "generic->Matriz4D_Transpose "), 1, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->mat4D_Transpose(&matin4d1);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = matin4d1;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->Matriz4D_Transpose  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestMatriz4DTransposeOf()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CMat4D tst,dst;
	const char *result;
	matin4d1.toZero();
	matin4d2.toZero();
	matReset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D TransposeOf()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->mat4D_TransposeOf(&matout4d,&matin4d1);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = matout4d;

		CProcClock::printClocks( va( "generic->Matriz4D_TransposeOf "), 1, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->mat4D_TransposeOf(&matout4d,&matin4d1);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = matout4d;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->Matriz4D_TransposeOf  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}



static void TestMatriz4DScale()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CMat4D tst[NUMTESTS],dst[NUMTESTS];
	const char *result;
	matin4d1.toZero();
	matin4d2.toZero();
	matReset4D();
	float value[NUMTESTS];
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D scale()..."<<endl;

	CRandom srnd( RANDOM_SEED );

	for ( int i = 0; i < NUMTESTS; i++ ) {
		value[i] = srnd.randomFloat();
		 dst[i].toZero();
		 tst[i].toZero();
	}

			for ( int i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			p_generic->mat4D_Scale(&matin4d1,value[i]);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
			tst[i] = matin4d1;
			}


		CProcClock::printClocks( va( "generic->Matriz4D_Scale "), NUMTESTS, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;



			for ( int i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			p_simd->mat4D_Scale(&matin4d1,value[i]);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
			dst[i] = matin4d1;
			}
		int i;
		for ( i = 0; i < NUMTESTS; i++ ) {
			if ( !dst[i].compare( tst[i], MATX_SIMD_EPSILON ) ) {
				break;
			}
		}
	result = ( i >= NUMTESTS ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->Matriz4D_Scale  %s",  result ), NUMTESTS, bestClocksSIMD, bestClocksGeneric );



}

static void TestMatriz4DScaleOf()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CMat4D tst[NUMTESTS],dst[NUMTESTS];
	const char *result;
	matin4d1.toZero();
	matin4d2.toZero();
	matReset4D();

	float value[NUMTESTS];

	CRandom srnd( RANDOM_SEED );

	for ( int i = 0; i < NUMTESTS; i++ ) {
		value[i] = srnd.randomFloat();
		 dst[i].toZero();
		 tst[i].toZero();
	}

	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D ScaleOf()..."<<endl;

			for ( int i = 0; i < NUMTESTS; i++ ) {

			StartRecordTimeLocal( start );
			p_generic->mat4D_ScaleOf(&matout4d,&matin4d1,value[i]);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
			tst[i] = matout4d;
			}
		CProcClock::printClocks( va( "generic->Matriz4D_ScaleOf "),NUMTESTS, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;



			for ( int i = 0; i < NUMTESTS; i++ ) {

			StartRecordTimeLocal( start );
			p_simd->mat4D_ScaleOf(&matout4d,&matin4d1,value[i]);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
			dst[i] = matout4d;
			}
		int i;
		for ( i = 0; i < NUMTESTS; i++ ) {
			if ( !dst[i].compare( tst[i], MATX_SIMD_EPSILON ) ) {
				break;
			}
		}
	result = ( i >= NUMTESTS ) ? "ok" : S_COLOR_RED"X";

	CProcClock::printClocks( va( "   simd->Matriz4D_ScaleOf  %s",  result ), NUMTESTS, bestClocksSIMD, bestClocksGeneric );



}



static void TestMatriz4DVectorMultiply()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec3D tst,dst;
	const char *result;
	matin4d1.toZero();
	matin4d2.toZero();
	matReset4D();
	reset3D();


	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D VectorMultiply()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->mat4D_VectorMultiply(&outVec,&matin4d1);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = outVec;

		CProcClock::printClocks( va( "generic->Matriz4D_VectorMultiply "), 1, bestClocksGeneric );
		matReset4D();
		reset3D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->mat4D_VectorMultiply(&outVec,&matin4d1);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = outVec;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->Matriz4D_VectorMultiply  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestMatriz4DVectorMultiplyOf()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec3D tst,dst;
	const char *result;
	matin4d1.toZero();
	matin4d2.toZero();
	matReset4D();
	reset3D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D VectorMultiplyOf()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->mat4D_VectorMultiplyOf(&outVec,&inVec1,&matin4d1);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = outVec;

		CProcClock::printClocks( va( "generic->Matriz4D_VectorMultiplyOf "), 1, bestClocksGeneric );
		matReset4D();
		reset3D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->mat4D_VectorMultiplyOf(&outVec,&inVec1,&matin4d1);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = outVec;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->Matriz4D_VectorMultiplyOf  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestMatriz4DVector4Multiply()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec4D tst,dst;
	const char *result;
	matin4d1.toZero();
	matin4d2.toZero();
	matReset4D();
	reset4D();


	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D VectorMultiply() CVec4D..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->mat4D_VectorMultiply(&outVec4d,&matin4d1);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = outVec4d;

		CProcClock::printClocks( va( "generic->Matriz4D_VectorMultiply (CVec4D) "), 1, bestClocksGeneric );
		matReset4D();
		reset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->mat4D_VectorMultiply(&outVec4d,&matin4d1);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = outVec4d;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->Matriz4D_VectorMultiply (CVec4D)  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestMatriz4DVector4MultiplyOf()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec4D tst,dst;
	const char *result;
	matin4d1.toZero();
	matin4d2.toZero();
	matReset4D();
	reset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D VectorMultiplyOf() (CVec4D)..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->mat4D_VectorMultiplyOf(&outVec4d,&inVec4d1,&matin4d1);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = outVec4d;

		CProcClock::printClocks( va( "generic->Matriz4D_VectorMultiplyOf (CVec4D)"), 1, bestClocksGeneric );
		matReset4D();
		reset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->mat4D_VectorMultiplyOf(&outVec4d,&inVec4d1,&matin4d1);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = outVec4d;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->Matriz4D_VectorMultiplyOf  (CVec4D) %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}


static void TestMatriz4DVectorToRotate()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CMat4D tst,dst;
	const char *result;
	matin4d1.toZero();
	matin4d2.toZero();
	matReset4D();
	float value[NUMTESTS];
	float value1[NUMTESTS];
	float value2[NUMTESTS];
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D ToRotate()..."<<endl;

	CRandom srnd( RANDOM_SEED );

	for ( int i = 0; i < NUMTESTS; i++ ) {
		value[i] = srnd.randomFloat();
		value1[i] = srnd.randomFloat();
		value2[i] = srnd.randomFloat();
		}

			for ( int i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			p_generic->mat4D_ToRotate(&matin4d1,value[i],value2[i],value2[i]);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
			}
		tst = matin4d1;

		CProcClock::printClocks( va( "generic->Matriz4D_ToRotate "), NUMTESTS, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;

		    dst.toZero();

			for ( int i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			p_simd->mat4D_ToRotate(&matin4d1,value[i],value2[i],value2[i]);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
			dst = matin4d1;
			}


		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->Matriz4D_ToRotate  %s",  result ), NUMTESTS, bestClocksSIMD, bestClocksGeneric );



}

static void TestMatriz4DVectorToRotateOf()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CMat4D tst,dst;
	const char *result;
	matin4d1.toZero();
	matin4d2.toZero();
	matReset4D();
	reset3D();
	CVec3D value[NUMTESTS];

	CRandom srnd( RANDOM_SEED );

	for ( int i = 0; i < NUMTESTS; i++ ) {
		value[i] = CVec3D(srnd.randomFloat(),srnd.randomFloat(),srnd.randomFloat());
		}

	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D ToRotateOf()..."<<endl;

			for ( int i = 0; i < NUMTESTS; i++ ) {

			StartRecordTimeLocal( start );
			p_generic->mat4D_ToRotateOf(&matin4d1,&value[i]);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
			tst = matin4d1;
			}
		CProcClock::printClocks( va( "generic->Matriz4D_ToRotateOf "),NUMTESTS, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;

		    dst.toZero();

			for ( int i = 0; i < NUMTESTS; i++ ) {

			StartRecordTimeLocal( start );
			p_simd->mat4D_ToRotateOf(&matin4d1,&value[i]);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
			dst = matin4d1;
			}
		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->Matriz4D_ToRotateOf  %s",  result ), NUMTESTS, bestClocksSIMD, bestClocksGeneric );



}



static void TestMatriz4DAlignedSum()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CMat4D tst,dst;
	const char *result;
	matin4d1Aligned.toZero();
	matin4d2Aligned.toZero();
	matReset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D AlignedSum()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->mat4D_AlignedSum(&matin4d1Aligned, &matin4d2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = matin4d1Aligned;

		CProcClock::printClocks( va( "generic->Matriz4D_AlignedSum "), 1, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->mat4D_AlignedSum(&matin4d1Aligned, &matin4d2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = matin4d1Aligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->Matriz4D_Sum  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestMatriz4DAlignedSumOf()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CMat4D tst,dst;
	const char *result;
	matin4d1Aligned.toZero();
	matin4d2Aligned.toZero();
	matReset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D AlignedSumOf()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->mat4D_AlignedSumOf(&matout4dAligned,&matin4d1Aligned, &matin4d2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = matin4d1Aligned;

		CProcClock::printClocks( va( "generic->Matriz4D_AlignedSumOf "), 1, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->mat4D_AlignedSumOf(&matout4dAligned,&matin4d1Aligned, &matin4d2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = matin4d1Aligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->Matriz4D_AlignedSumOf  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestMatriz4DAlignedDiff()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CMat4D tst,dst;
	const char *result;
	matin4d1Aligned.toZero();
	matin4d2Aligned.toZero();
	matReset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D AlignedDiff()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->mat4D_AlignedDiff(&matin4d1Aligned, &matin4d2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = matin4d1Aligned;

		CProcClock::printClocks( va( "generic->Matriz4D_AlignedDiff "), 1, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->mat4D_AlignedDiff(&matin4d1Aligned, &matin4d2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = matin4d1Aligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->Matriz4D_AlignedDiff  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestMatriz4DAlignedDiffOf()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CMat4D tst,dst;
	const char *result;
	matin4d1Aligned.toZero();
	matin4d2Aligned.toZero();
	matReset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D AlignedDiffOf()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->mat4D_AlignedDiffOf(&matout4dAligned,&matin4d1Aligned, &matin4d2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = matout4dAligned;

		CProcClock::printClocks( va( "generic->Matriz4D_AlignedDiffOf "), 1, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->mat4D_AlignedDiffOf(&matout4dAligned,&matin4d1Aligned, &matin4d2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = matout4dAligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->Matriz4D_AlignedDiffOf  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestMatriz4DAlignedMultiply()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CMat4D tst,dst;
	const char *result;
	matin4d1Aligned.toZero();
	matin4d2Aligned.toZero();
	matReset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D AlignedMultiply()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->mat4D_AlignedMultiply(&matin4d1Aligned, &matin4d2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = matin4d1Aligned;

		CProcClock::printClocks( va( "generic->Matriz4D_AlignedMultiply "), 1, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->mat4D_AlignedMultiply(&matin4d1Aligned, &matin4d2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = matin4d1Aligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->Matriz4D_AlignedMultiply  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestMatriz4DAlignedMultiplyOf()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CMat4D tst,dst;
	const char *result;
	matin4d1Aligned.toZero();
	matin4d2Aligned.toZero();
	matReset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D AlignedMultiplyOf()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->mat4D_AlignedMultiplyOf(&matout4dAligned,&matin4d1Aligned, &matin4d2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = matout4dAligned;

		CProcClock::printClocks( va( "generic->Matriz4D_AlignedMultiplyOf "), 1, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->mat4D_AlignedMultiplyOf(&matout4dAligned,&matin4d1Aligned, &matin4d2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = matout4dAligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->Matriz4D_AlignedMultiplyOf  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}



static void TestMatriz4DAlignedTranspose()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CMat4D tst,dst;
	const char *result;
	matin4d1Aligned.toZero();
	matin4d2Aligned.toZero();
	matReset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D AlignedTranspose()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->mat4D_AlignedTranspose(&matin4d1Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = matin4d1Aligned;

		CProcClock::printClocks( va( "generic->Matriz4D_AlignedTranspose "), 1, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->mat4D_AlignedTranspose(&matin4d1Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = matin4d1Aligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->Matriz4D_AlignedTranspose  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestMatriz4DAlignedTransposeOf()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CMat4D tst,dst;
	const char *result;
	matin4d1Aligned.toZero();
	matin4d2Aligned.toZero();
	matReset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D AlignedTransposeOf()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->mat4D_AlignedTransposeOf(&matout4dAligned,&matin4d1Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = matout4dAligned;

		CProcClock::printClocks( va( "generic->Matriz4D_AlignedTransposeOf "), 1, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->mat4D_AlignedTransposeOf(&matout4dAligned,&matin4d1Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = matout4dAligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->Matriz4D_AlignedTransposeOf  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}



static void TestMatriz4DAlignedScale()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CMat4D tst[NUMTESTS],dst[NUMTESTS];
	const char *result;
	matin4d1Aligned.toZero();
	matin4d2Aligned.toZero();
	matReset4D();
	float value[NUMTESTS];
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D AlignedScale()..."<<endl;

	CRandom srnd( RANDOM_SEED );

	for ( int i = 0; i < NUMTESTS; i++ ) {
		value[i] = srnd.randomFloat();
		dst[i].toZero();
		tst[i].toZero();
		}

			for ( int i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			p_generic->mat4D_AlignedScale(&matin4d1Aligned,value[i]);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
			tst[i] = matin4d1Aligned;
			}

		CProcClock::printClocks( va( "generic->Matriz4D_AlignedScale "), NUMTESTS, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;



			for ( int i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			p_simd->mat4D_AlignedScale(&matin4d1Aligned,value[i]);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
			dst[i] = matin4d1Aligned;
			}
				int i;
		for ( i = 0; i < NUMTESTS; i++ ) {
			if ( !dst[i].compare( tst[i], MATX_SIMD_EPSILON ) ) {
				break;
			}
		}
		result = ( i >= NUMTESTS ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->Matriz4D_AlignedScale  %s",  result ), NUMTESTS, bestClocksSIMD, bestClocksGeneric );



}

static void TestMatriz4DAlignedScaleOf()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CMat4D tst[NUMTESTS],dst[NUMTESTS];
	const char *result;
	matin4d1.toZero();
	matin4d2.toZero();
	matReset4D();

	float value[NUMTESTS];

	CRandom srnd( RANDOM_SEED );

	for ( int i = 0; i < NUMTESTS; i++ ) {
		value[i] = srnd.randomFloat();
		 dst[i].toZero();
		 tst[i].toZero();
	}

	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D AlignedScaleOf()..."<<endl;

			for ( int i = 0; i < NUMTESTS; i++ ) {
			matReset4D();
			StartRecordTimeLocal( start );
			p_generic->mat4D_AlignedScaleOf(&matout4dAligned,&matin4d1Aligned,value[i]);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
			tst[i] = matout4dAligned;
			}
		CProcClock::printClocks( va( "generic->Matriz4D_AlignedScaleOf "),NUMTESTS, bestClocksGeneric );

		bestClocksSIMD = 0;



			for ( int i = 0; i < NUMTESTS; i++ ) {
			matReset4D();
			StartRecordTimeLocal( start );
			p_simd->mat4D_AlignedScaleOf(&matout4dAligned,&matin4d1Aligned,value[i]);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
			dst[i] = matout4dAligned;
			}
		int i;
		for ( i = 0; i < NUMTESTS; i++ ) {
			if ( !dst[i].compare( tst[i], MATX_SIMD_EPSILON ) ) {
				break;
			}
		}
	result = ( i >= NUMTESTS ) ? "ok" : S_COLOR_RED"X";

	CProcClock::printClocks( va( "   simd->Matriz4D_AlignedScaleOf  %s",  result ), NUMTESTS, bestClocksSIMD, bestClocksGeneric );



}


static void TestMatriz4DAlignedVectorMultiply()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec3D tst,dst;
	const char *result;
	matin4d1.toZero();
	matin4d2Aligned.toZero();
	matReset4D();
	reset3D();


	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D AlignedVectorMultiply()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->mat4D_AlignedVectorMultiply(&outVecAligned,&matin4d1Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = outVecAligned;

		CProcClock::printClocks( va( "generic->Matriz4D_AlignedVectorMultiply "), 1, bestClocksGeneric );
		matReset4D();
		reset3D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->mat4D_AlignedVectorMultiply(&outVecAligned,&matin4d1Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = outVecAligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->Matriz4D_AlignedVectorMultiply  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}



static void TestMatriz4DAlignedVectorMultiplyOf()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec3D tst,dst;
	const char *result;
	matin4d1Aligned.toZero();
	matin4d2Aligned.toZero();
	matReset4D();
	reset3D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D AlignedVectorMultiplyOf()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->mat4D_AlignedVectorMultiplyOf(&outVecAligned,&inVec1Aligned,&matin4d1Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = outVecAligned;

		CProcClock::printClocks( va( "generic->Matriz4D_AlignedVectorMultiplyOf "), 1, bestClocksGeneric );
		matReset4D();
		reset3D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->mat4D_AlignedVectorMultiplyOf(&outVecAligned,&inVec1Aligned,&matin4d1Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = outVecAligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->Matriz4D_AlignedVectorMultiplyOf  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestMatriz4DAlignedVector4Multiply()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec4D tst,dst;
	const char *result;
	matin4d1Aligned.toZero();
	matin4d2Aligned.toZero();
	matReset4D();
	reset4D();


	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D AlignedVectorMultiply() CVec4D..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->mat4D_AlignedVectorMultiply(&outVec4dAligned,&matin4d1Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = outVec4dAligned;

		CProcClock::printClocks( va( "generic->Matriz4D_AlignedVectorMultiply (CVec4D) "), 1, bestClocksGeneric );
		matReset4D();
		reset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->mat4D_AlignedVectorMultiply(&outVec4dAligned,&matin4d1Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = outVec4dAligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->Matriz4D_AlignedVectorMultiply (CVec4D)  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestMatriz4DAlignedVector4MultiplyOf()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec4D tst,dst;
	const char *result;
	matin4d1Aligned.toZero();
	matin4d2Aligned.toZero();
	matReset4D();
	reset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D AlignedVectorMultiplyOf() (CVec4D)..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->mat4D_AlignedVectorMultiplyOf(&outVec4dAligned,&inVec4d1Aligned,&matin4d1Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = outVec4dAligned;

		CProcClock::printClocks( va( "generic->Matriz4D_AlignedVectorMultiplyOf (CVec4D)"), 1, bestClocksGeneric );
		matReset4D();
		reset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->mat4D_AlignedVectorMultiplyOf(&outVec4dAligned,&inVec4d1Aligned,&matin4d1Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = outVec4dAligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->Matriz4D_AlignedVectorMultiplyOf  (CVec4D) %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}


static void TestMatriz4DAlignedVectorToRotate()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CMat4D tst,dst;
	const char *result;
	matin4d1Aligned.toZero();
	matin4d2Aligned.toZero();
	matReset4D();
	float value[NUMTESTS];
	float value1[NUMTESTS];
	float value2[NUMTESTS];
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D AlignedToRotate()..."<<endl;

	CRandom srnd( RANDOM_SEED );

	for ( int i = 0; i < NUMTESTS; i++ ) {
		value[i] = srnd.randomFloat();
		value1[i] = srnd.randomFloat();
		value2[i] = srnd.randomFloat();
		}

			for ( int i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			p_generic->mat4D_AlignedToRotate(&matin4d1Aligned,value[i],value2[i],value2[i]);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
			}
		tst = matin4d1Aligned;

		CProcClock::printClocks( va( "generic->Matriz4D_AlignedToRotate "), NUMTESTS, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;

		    dst.toZero();

			for ( int i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			p_simd->mat4D_AlignedToRotate(&matin4d1Aligned,value[i],value2[i],value2[i]);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
			dst = matin4d1Aligned;
			}


		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->Matriz4D_AlignedToRotate  %s",  result ), NUMTESTS, bestClocksSIMD, bestClocksGeneric );



}

static void TestMatriz4DAlignedVectorToRotateOf()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CMat4D tst,dst;
	const char *result;
	matin4d1Aligned.toZero();
	matin4d2Aligned.toZero();
	matReset4D();
	reset3D();
	CVec3D value[NUMTESTS];

	CRandom srnd( RANDOM_SEED );

	for ( int i = 0; i < NUMTESTS; i++ ) {
		value[i] = CVec3D(srnd.randomFloat(),srnd.randomFloat(),srnd.randomFloat());
		}

	Debug::debug(Debug::math,__FUNCTION__) << "Testing Matriz4D AlignedToRotateOf()..."<<endl;

			for ( int i = 0; i < NUMTESTS; i++ ) {

			StartRecordTimeLocal( start );
			p_generic->mat4D_AlignedToRotateOf(&matin4d1Aligned,&value[i]);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
			tst = matin4d1Aligned;
			}
		CProcClock::printClocks( va( "generic->Matriz4D_AlignedToRotateOf "),NUMTESTS, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;

		    dst.toZero();

			for ( int i = 0; i < NUMTESTS; i++ ) {

			StartRecordTimeLocal( start );
			p_simd->mat4D_AlignedToRotateOf(&matin4d1Aligned,&value[i]);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
			dst = matin4d1Aligned;
			}
		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->Matriz4D_AlignedToRotateOf  %s",  result ), NUMTESTS, bestClocksSIMD, bestClocksGeneric );



}





void CSIMD::TestMat4D_f( ){

#ifdef _WIN32
	SetThreadPriority( GetCurrentThread(), THREAD_PRIORITY_TIME_CRITICAL );
#endif /* _WIN32 */

	p_simd = processor;
	p_generic = generic;



	//idLib::common->SetRefreshOnPrint( true );

	Debug::debug(Debug::math,__FUNCTION__) << "using %s for SIMD processing: "<< p_simd->getName() <<endl;

	CProcClock::getBaseClocks();
	//testes=============================================

	TestMatriz4DSum();
	TestMatriz4DSumOf();

	TestMatriz4DDiff();
	TestMatriz4DDiffOf();
	TestMatriz4DMultiply();
	TestMatriz4DMultiplyOf();
	TestMatriz4DTranspose();
	TestMatriz4DTransposeOf();
	TestMatriz4DScale();
	TestMatriz4DScaleOf();
	TestMatriz4DVectorMultiply();
	TestMatriz4DVectorMultiplyOf();
	TestMatriz4DVector4Multiply();
	TestMatriz4DVector4MultiplyOf();
	TestMatriz4DVectorToRotate();
	TestMatriz4DVectorToRotateOf();



	TestMatriz4DAlignedSum();
	TestMatriz4DAlignedSumOf();
	TestMatriz4DAlignedDiff();
	TestMatriz4DAlignedDiffOf();
	TestMatriz4DAlignedScale();
	TestMatriz4DAlignedScaleOf();
	TestMatriz4DAlignedMultiply();
	TestMatriz4DAlignedMultiplyOf();
	TestMatriz4DAlignedTranspose();
	TestMatriz4DAlignedTransposeOf();
	TestMatriz4DAlignedVectorMultiply();
	TestMatriz4DAlignedVectorMultiplyOf();
	TestMatriz4DAlignedVector4Multiply();
	TestMatriz4DAlignedVector4MultiplyOf();
	TestMatriz4DAlignedVectorToRotate();
	TestMatriz4DAlignedVectorToRotateOf();

	// fim testes==========================================
	if ( p_simd != processor ) {
		delete p_simd;
	}
	p_simd = NULL;
	p_generic = NULL;

#ifdef _WIN32
	SetThreadPriority( GetCurrentThread(), THREAD_PRIORITY_NORMAL );
#endif /* _WIN32 */

}









} //end math
} //end SMF
