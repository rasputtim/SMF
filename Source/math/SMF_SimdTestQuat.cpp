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
#include "math/SMF_Quaternion.h"
#include "math/SMF_EulerAngles.h"
#include "sys/SMF_System.h"
#include "util/SMF_Debug.h"
#include "util/SMF_StringUtils.h"

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
extern CSIMDProcessor *	generic;			// pointer to generic SIMD implementation
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




//===================CQuaternion====================================



CQuaternion matout4d, matin4d1, matin4d2;
ALIGNTO16  CQuaternion matout4dAligned;
ALIGNTO16  CQuaternion matin4d1Aligned;
ALIGNTO16  CQuaternion matin4d2Aligned;
static void matReset4D()
{
	matout4d.set(0,0,0,0);
	matout4dAligned.set(0,0,0,0);
	CEulerAngles ang1(90,45,0);
	CEulerAngles ang2(0,90,45);

	matin4d1 = ang1.toQuat();
	matin4d2 = ang2.toQuat();

	matin4d1Aligned = matin4d1;
	matin4d2Aligned = matin4d2;

}



static void Test_quaternion_Normalize()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CQuaternion tst[NUMTESTS],dst[NUMTESTS];
	const char *result;
	matReset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing  quaternion_Normalize..."<<endl;

	CRandom srnd( RANDOM_SEED );

	for ( int i = 0; i < NUMTESTS; i++ ) {
		dst[i].set(0,0,0,0);
		tst[i].set(0,0,0,0);
	}

			for ( int i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			p_generic-> quaternion_Normalize(&matin4d1);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
			tst[i] = matin4d1;
			}


		CProcClock::printClocks( va( "generic-> quaternion_Normalize "), NUMTESTS, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;



			for ( int i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			p_simd-> quaternion_Normalize(&matin4d1);
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
	CProcClock::printClocks( va( "   simd-> quaternion_Normalize  %s",  result ), NUMTESTS, bestClocksSIMD, bestClocksGeneric );



}


static void Test_quaternion_NormalizeOf()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CQuaternion tst[NUMTESTS],dst[NUMTESTS];
	const char *result;
	matReset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing  quaternion_NormalizeOf..."<<endl;

	CRandom srnd( RANDOM_SEED );

	for ( int i = 0; i < NUMTESTS; i++ ) {
		dst[i].set(0,0,0,0);
		tst[i].set(0,0,0,0);
	}

			for ( int i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			p_generic-> quaternion_NormalizeOf(&matout4d,&matin4d1);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
			tst[i] = matout4d;
			}


		CProcClock::printClocks( va( "generic-> quaternion_NormalizeOf "), NUMTESTS, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;



			for ( int i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			p_simd-> quaternion_NormalizeOf(&matout4d,&matin4d1);
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
	CProcClock::printClocks( va( "   simd-> quaternion_NormalizeOf  %s",  result ), NUMTESTS, bestClocksSIMD, bestClocksGeneric );



}

static void Test_quaternion_Multiply()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CQuaternion tst[NUMTESTS],dst[NUMTESTS];
	const char *result;
	matReset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing  quaternion_Multiply..."<<endl;

	CRandom srnd( RANDOM_SEED );

	for ( int i = 0; i < NUMTESTS; i++ ) {
		dst[i].set(0,0,0,0);
		tst[i].set(0,0,0,0);
	}

			for ( int i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			p_generic-> quaternion_Multiply(&matin4d1,&matin4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
			tst[i] = matin4d1;
			}


		CProcClock::printClocks( va( "generic-> quaternion_Multiply "), NUMTESTS, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;



			for ( int i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			p_simd-> quaternion_Multiply(&matin4d1,&matin4d2);
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
	CProcClock::printClocks( va( "   simd-> quaternion_Multiply  %s",  result ), NUMTESTS, bestClocksSIMD, bestClocksGeneric );



}


static void Test_quaternion_MultiplyOf()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CQuaternion tst[NUMTESTS],dst[NUMTESTS];
	const char *result;
	matReset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing  quaternion_MultiplyOf..."<<endl;

	CRandom srnd( RANDOM_SEED );

	for ( int i = 0; i < NUMTESTS; i++ ) {
		dst[i].set(0,0,0,0);
		tst[i].set(0,0,0,0);
	}

			for ( int i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			p_generic-> quaternion_MultiplyOf(&matout4d,&matin4d1,&matin4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
			tst[i] = matout4d;
			}


		CProcClock::printClocks( va( "generic-> quaternion_MultiplyOf "), NUMTESTS, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;



			for ( int i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			p_simd-> quaternion_MultiplyOf(&matout4d,&matin4d1,&matin4d2);
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
	CProcClock::printClocks( va( "   simd-> quaternion_MultiplyOf  %s",  result ), NUMTESTS, bestClocksSIMD, bestClocksGeneric );



}


void CSIMD::TestQuat4D_f(){

#ifdef _WIN32
	SetThreadPriority( GetCurrentThread(), THREAD_PRIORITY_TIME_CRITICAL );
#endif /* _WIN32 */

	p_simd = processor;
	p_generic = generic;



	//idLib::common->SetRefreshOnPrint( true );

	Debug::debug(Debug::math,__FUNCTION__) << "using %s for SIMD processing: "<< p_simd->getName() <<endl;

	CProcClock::getBaseClocks();
	//testes=============================================

	Test_quaternion_Normalize();
	Test_quaternion_NormalizeOf();

	Test_quaternion_Multiply();
	Test_quaternion_MultiplyOf();


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
