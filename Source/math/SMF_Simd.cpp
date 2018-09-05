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
#include "math/SMF_SimdSSE41.h"
#include "math/SMF_SimdSSE42.h"
#include "math/SMF_SimdAVX.h"
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
using namespace Util;
namespace MATH{
using namespace GEO;

const sf_m128 sseMaskXYZ = _mm_set_ps(0.f, andMaskOneF, andMaskOneF, andMaskOneF);

const sf_m128 sseSignMask = _mm_set1_ps(-0.f); // -0.f = 1 << 31
const sf_m128 sseSignMask3 = _mm_set_ps(0.f, -0.f, -0.f, -0.f); // -0.f = 1 << 31
const sf_m128 sseSignMask4 = _mm_set_ps(1.0f, 0.f, 0.f, 0.f); // -0.f = 1 << 31
const __m256 sseSignMask256 = _mm256_set1_ps(-0.f);
CSIMDProcessor	*	processor = NULL;			// pointer to SIMD processor
CSIMDProcessor *	generic = NULL;				// pointer to generic SIMD implementation
CSIMDProcessor *	SIMDProcessor = NULL;

#define TIME_TYPE unsigned int
#define	MATX_SIMD_EPSILON 1e-5f

sf_m128 FLOAT_TO_M128(float f)
{
	return _mm_load_ss(&f);
}
float M128_TO_FLOAT(sf_m128 sse)
{
	float ret;
	_mm_store_ss(&ret, sse);
	return ret;
}
sf_m128 SET_PS(float um,float dois, float tres, float quatro){
return _mm_set_ps(um,dois,tres,quatro);
}
/*
================
CSIMD::init
================
*/
void CSIMD::init() {
	generic = new CSIMD_Generic;
	generic->cpuid = CPUID_GENERIC;
	processor = NULL;
	SIMDProcessor = generic;
}
bool CSIMD::isInitialized(){
	if(generic == NULL) return false;
	else return true;
}

bool CSIMD::globalUseSIMD(){
	return true;//CGlobalConfiguration::useSIMDprocessor.getBool();
}
void  CSIMD::initHeap(){
	if(!mem_isInitialized()) mem_Init();
	if(!MATH::CMath::isInitialized()) MATH::CMath::init();
	if(!MATH::CSIMD::isInitialized()) MATH::CSIMD::init();
	if(!MATH::CSIMD::isProcessorInitialized()) MATH::CSIMD::initProcessor( "sgfInternal", false );

}
/*
============
CSIMD::InitProcessor
============
*/

CSIMDProcessor *CSIMD::getProcessor(){
	if (processor==NULL) initProcessor("",false);
	return processor;
}

CSIMDProcessor *CSIMD::getGenProcessor(){
	if (generic==NULL) init();
	return generic;
}

void CSIMD::initProcessor( const char *module, bool forceGeneric ) {
	cpuid_t cpuid;
	CSIMDProcessor *newProcessor;

	cpuid = System::getProcessorId();

	if ( forceGeneric ) {

		newProcessor = generic;
		processor = generic;

	} else {

		if ( !processor ) {
			if ( ( cpuid & CPUID_ALTIVEC ) ) {
				//processor = new CSIMD_AltiVec;
				processor = generic;
			}else if ( ( cpuid & CPUID_MMX ) && ( cpuid & CPUID_SSE ) && ( cpuid & CPUID_SSE2 ) && ( cpuid & CPUID_SSE3 ) && ( cpuid & CPUID_SSE41 )&& ( cpuid & CPUID_SSE42 )&& ( cpuid & CPUID_AVX )) {
				processor = new CSIMD_AVX;//CSIMD_SSE42;
			}else if ( ( cpuid & CPUID_MMX ) && ( cpuid & CPUID_SSE ) && ( cpuid & CPUID_SSE2 ) && ( cpuid & CPUID_SSE3 ) && ( cpuid & CPUID_SSE41 )&& ( cpuid & CPUID_SSE42 )) {
				processor = new CSIMD_SSE42;//CSIMD_SSE42;
			}else if ( ( cpuid & CPUID_MMX ) && ( cpuid & CPUID_SSE ) && ( cpuid & CPUID_SSE2 ) && ( cpuid & CPUID_SSE3 ) && ( cpuid & CPUID_SSE41 )) {
				processor = new CSIMD_SSE41;//CCPUID_SSE41;
			} else if ( ( cpuid & CPUID_MMX ) && ( cpuid & CPUID_SSE ) && ( cpuid & CPUID_SSE2 ) && ( cpuid & CPUID_SSE3 ) ) {
				processor = new CSIMD_SSE3;
			} else if ( ( cpuid & CPUID_MMX ) && ( cpuid & CPUID_SSE ) && ( cpuid & CPUID_SSE2 ) ) {
				processor = new CSIMD_SSE2;
			} else if ( ( cpuid & CPUID_MMX ) && ( cpuid & CPUID_SSE ) ) {
				processor = new CSIMD_SSE;
			#if !defined (__GNUC__)
			} else if ( ( cpuid & CPUID_MMX ) && ( cpuid & CPUID_3DNOW ) ) {
				processor = new CSIMD_3DNow;
            #endif
			} else if ( ( cpuid & CPUID_MMX ) ) {
				processor = new CSIMD_MMX;
			} else {
				processor = generic;
			}
			processor->cpuid = cpuid;
		}

		newProcessor = processor;
	}

	if ( newProcessor != SIMDProcessor ) {
		SIMDProcessor = newProcessor;
		Debug::debug(Debug::math,__FUNCTION__) <<  module << " using "  << SIMDProcessor->getName()<<" for SIMD processing"<< endl;
	}
#if defined (_MSC_VER) //mingw cause segmentation fault when ftz or daz are enables toguether with HIPREC

	if ( cpuid & CPUID_FTZ ) {
		System::CFPU::fpu_SetFTZ( true );
		Debug::debug(Debug::math,__FUNCTION__) << "enabled Flush-To-zero mode" <<endl;
	}

	if ( cpuid & CPUID_DAZ ) {
		System::CFPU::fpu_SetDAZ( true );
		Debug::debug(Debug::math,__FUNCTION__) << "enabled Denormals-Are-zero mode\n" <<endl;
	}
#endif
}

bool		CSIMD::isProcessorInitialized(){
	if(processor==NULL) return false;
	else return true;
}
/*
================
CSIMD::shutdown
================
*/
void CSIMD::shutdown() {
	mem_Shutdown();
	if ( processor != generic ) {
		delete processor;
	}
	delete generic;
	generic = NULL;
	processor = NULL;
	SIMDProcessor = NULL;

}


//===============================================================
//
// Test code
//
//===============================================================

#define COUNT		1024		// data count
#define NUMTESTS	2048		// number of tests

#define RANDOM_SEED		1013904223L	//((int)idLib::sys->GetClockTicks())

CSIMDProcessor *p_simd;
CSIMDProcessor *p_generic;
unsigned baseClocks = 0;

unsigned saved_ebx = 0;
unsigned start_ClockCount = 0;
unsigned end_ClockCount = 0;
double ticksPerNanosecond;

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



#define GetBest( start, end, best )			\
	if ( !best || end - start < best ) {	\
		best = end - start;					\
	}

#endif
/*
============
CProcClock::printClocks
============
*/
void CProcClock::printClocks( const char *string, int dataCount, int clocks, int otherClocks ) {
	int i;

	Debug::debug(Debug::math,__FUNCTION__) << string <<endl;
	//for ( i = CMyString::lengthWithoutColors(string); i < 48; i++ ) {
	//	Debug::debug(Debug::math,__FUNCTION__) <<" "<<endl;
	//}
	clocks -= baseClocks;
	if ( otherClocks && clocks ) {
		otherClocks -= baseClocks;
		int p = (int) ( (float) ( otherClocks - clocks ) * 100.0f / (float) otherClocks );
		Debug::debug(Debug::math,__FUNCTION__) << "count = "<< dataCount<<" , clocks =  "<< clocks <<", "<<  p<< "% melhor"<<endl;
	} else {
		Debug::debug(Debug::math,__FUNCTION__) << "count =  "<<  dataCount<< ", clocks = "<< clocks <<endl;
	}
}

/*
============
GetBaseClocks
============
*/
void CProcClock::getBaseClocks() {
	int i, start, end, bestClocks;

	bestClocks = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
	}
	baseClocks = bestClocks;
}

/*
============
TestAdd
============
*/
void TestAdd() {
	int i;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( float fdst0[COUNT] );
	ALIGN16( float fdst1[COUNT] );
	ALIGN16( float fsrc0[COUNT] );
	ALIGN16( float fsrc1[COUNT] );
	const char *result;

	CRandom srnd( RANDOM_SEED );

	for ( i = 0; i < COUNT; i++ ) {
		fsrc0[i] = srnd.randomFloat() * 10.0f;
		fsrc1[i] = srnd.randomFloat() * 10.0f;
	}

	Debug::debug(Debug::math,__FUNCTION__) <<"====================================" <<endl;

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->add( fdst0, 4.0f, fsrc1, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->add( float + float[] )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->add( fdst1, 4.0f, fsrc1, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( CMath::fabs( fdst0[i] - fdst1[i] ) > 1e-5f ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->add( float + float[] ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->add( fdst0, fsrc0, fsrc1, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->add( float[] + float[] )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->add( fdst1, fsrc0, fsrc1, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( CMath::fabs( fdst0[i] - fdst1[i] ) > 1e-5f ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->add( float[] + float[] ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );
}

/*
============
TestSub
============
*/
void TestSub() {
	int i;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( float fdst0[COUNT] );
	ALIGN16( float fdst1[COUNT] );
	ALIGN16( float fsrc0[COUNT] );
	ALIGN16( float fsrc1[COUNT] );
	const char *result;

	CRandom srnd( RANDOM_SEED );

	for ( i = 0; i < COUNT; i++ ) {
		fsrc0[i] = srnd.randomFloat() * 10.0f;
		fsrc1[i] = srnd.randomFloat() * 10.0f;
	}

	Debug::debug(Debug::math,__FUNCTION__) <<"====================================\n" <<endl;

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->sub( fdst0, 4.0f, fsrc1, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->sub( float + float[] )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->sub( fdst1, 4.0f, fsrc1, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( CMath::fabs( fdst0[i] - fdst1[i] ) > 1e-5f ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->sub( float + float[] ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->sub( fdst0, fsrc0, fsrc1, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->sub( float[] + float[] )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->sub( fdst1, fsrc0, fsrc1, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( CMath::fabs( fdst0[i] - fdst1[i] ) > 1e-5f ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->sub( float[] + float[] ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );
}

/*
============
TestMul
============
*/
void TestMul() {
	int i;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( float fdst0[COUNT] );
	ALIGN16( float fdst1[COUNT] );
	ALIGN16( float fsrc0[COUNT] );
	ALIGN16( float fsrc1[COUNT] );
	const char *result;

	CRandom srnd( RANDOM_SEED );

	for ( i = 0; i < COUNT; i++ ) {
		fsrc0[i] = srnd.randomFloat() * 10.0f;
		fsrc1[i] = srnd.randomFloat() * 10.0f;
	}

	Debug::debug(Debug::math,__FUNCTION__) <<"====================================" <<endl;

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->mul( fdst0, 4.0f, fsrc1, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->mul( float * float[] )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->mul( fdst1, 4.0f, fsrc1, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( CMath::fabs( fdst0[i] - fdst1[i] ) > 1e-5f ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->mul( float * float[] ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );


	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->mul( fdst0, fsrc0, fsrc1, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->mul( float[] * float[] )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->mul( fdst1, fsrc0, fsrc1, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( CMath::fabs( fdst0[i] - fdst1[i] ) > 1e-5f ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->mul( float[] * float[] ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );
}

/*
============
TestDiv
============
*/
void TestDiv() {
	int i;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( float fdst0[COUNT] );
	ALIGN16( float fdst1[COUNT] );
	ALIGN16( float fsrc0[COUNT] );
	ALIGN16( float fsrc1[COUNT] );
	const char *result;

	CRandom srnd( RANDOM_SEED );

	for ( i = 0; i < COUNT; i++ ) {
		fsrc0[i] = srnd.randomFloat() * 10.0f;
		do {
			fsrc1[i] = srnd.randomFloat() * 10.0f;
		} while( CMath::fabs( fsrc1[i] ) < 0.1f );
	}

	Debug::debug(Debug::math,__FUNCTION__) <<"===================================="<<endl;


	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->div( fdst0, 4.0f, fsrc1, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->div( float * float[] )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->div( fdst1, 4.0f, fsrc1, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( CMath::fabs( fdst0[i] - fdst1[i] ) > 1e-5f ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->div( float * float[] ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );


	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->div( fdst0, fsrc0, fsrc1, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->div( float[] * float[] )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->div( fdst1, fsrc0, fsrc1, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( CMath::fabs( fdst0[i] - fdst1[i] ) > 1e-3f ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->div( float[] * float[] ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );
}

/*
============
TestMulAdd
============
*/
void TestMulAdd() {
	int i, j;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( float fdst0[COUNT] );
	ALIGN16( float fdst1[COUNT] );
	ALIGN16( float fsrc0[COUNT] );
	const char *result;

	CRandom srnd( RANDOM_SEED );

	for ( i = 0; i < COUNT; i++ ) {
		fsrc0[i] = srnd.randomFloat() * 10.0f;
	}

	Debug::debug(Debug::math,__FUNCTION__) <<"====================================" <<endl;

	for ( j = 0; j < 50 && j < COUNT; j++ ) {

		bestClocksGeneric = 0;
		for ( i = 0; i < NUMTESTS; i++ ) {
			for ( int k = 0; k < COUNT; k++ ) {
				fdst0[k] = k;
			}
			StartRecordTimeLocal( start );
			p_generic->mulAdd( fdst0, 0.123f, fsrc0, j );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
		}
		CProcClock::printClocks( va( "generic->mulAdd( float * float[%2d] )", j ), 1, bestClocksGeneric );

		bestClocksSIMD = 0;
		for ( i = 0; i < NUMTESTS; i++ ) {
			for ( int k = 0; k < COUNT; k++ ) {
				fdst1[k] = k;
			}
			StartRecordTimeLocal( start );
			p_simd->mulAdd( fdst1, 0.123f, fsrc0, j );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		}

		for ( i = 0; i < COUNT; i++ ) {
			if ( CMath::fabs( fdst0[i] - fdst1[i] ) > 1e-5f ) {
				break;
			}
		}
		result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->mulAdd( float * float[%2d] ) %s", j, result ), 1, bestClocksSIMD, bestClocksGeneric );
	}
}

/*
============
TestMulSub
============
*/
void TestMulSub() {
	int i, j;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( float fdst0[COUNT] );
	ALIGN16( float fdst1[COUNT] );
	ALIGN16( float fsrc0[COUNT] );
	const char *result;

	CRandom srnd( RANDOM_SEED );

	for ( i = 0; i < COUNT; i++ ) {
		fsrc0[i] = srnd.randomFloat() * 10.0f;
	}

	Debug::debug(Debug::math,__FUNCTION__) <<"====================================" <<endl;

	for ( j = 0; j < 50 && j < COUNT; j++ ) {

		bestClocksGeneric = 0;
		for ( i = 0; i < NUMTESTS; i++ ) {
			for ( int k = 0; k < COUNT; k++ ) {
				fdst0[k] = k;
			}
			StartRecordTimeLocal( start );
			p_generic->mulSub( fdst0, 0.123f, fsrc0, j );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
		}
		CProcClock::printClocks( va( "generic->mulSub( float * float[%2d] )", j ), 1, bestClocksGeneric );

		bestClocksSIMD = 0;
		for ( i = 0; i < NUMTESTS; i++ ) {
			for ( int k = 0; k < COUNT; k++ ) {
				fdst1[k] = k;
			}
			StartRecordTimeLocal( start );
			p_simd->mulSub( fdst1, 0.123f, fsrc0, j );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		}

		for ( i = 0; i < COUNT; i++ ) {
			if ( CMath::fabs( fdst0[i] - fdst1[i] ) > 1e-5f ) {
				break;
			}
		}
		result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->mulSub( float * float[%2d] ) %s", j, result ), 1, bestClocksSIMD, bestClocksGeneric );
	}
}

/*
============
TestDot
============
*/
void TestDot() {
	int i, j;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( float fdst0[COUNT] );
	ALIGN16( float fdst1[COUNT] );
	ALIGN16( float fsrc0[COUNT] );
	ALIGN16( float fsrc1[COUNT] );
	ALIGN16( CVec3D v3src0[COUNT] );
	ALIGN16( CVec3D v3src1[COUNT] );
	ALIGN16( CVec3D v3constant ) ( 1.0f, 2.0f, 3.0f );
	ALIGN16( GEO::CPlane v4src0[COUNT] );
	ALIGN16( GEO::CPlane v4constant ) (1.0f, 2.0f, 3.0f, 4.0f);
	ALIGN16( SMF:: CVertex drawVerts[COUNT] );
	const char *result;

	CRandom srnd( RANDOM_SEED );

	for ( i = 0; i < COUNT; i++ ) {
		fsrc0[i] = srnd.randomFloat() * 10.0f;
		fsrc1[i] = srnd.randomFloat() * 10.0f;
		v3src0[i][0] = srnd.randomFloat() * 10.0f;
		v3src0[i][1] = srnd.randomFloat() * 10.0f;
		v3src0[i][2] = srnd.randomFloat() * 10.0f;
		v3src1[i][0] = srnd.randomFloat() * 10.0f;
		v3src1[i][1] = srnd.randomFloat() * 10.0f;
		v3src1[i][2] = srnd.randomFloat() * 10.0f;
		v4src0[i] = v3src0[i];
		v4src0[i][3] = srnd.randomFloat() * 10.0f;
		drawVerts[i].xyz = v3src0[i];
	}

	Debug::debug(Debug::math,__FUNCTION__) <<"====================================" <<endl;


	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->dot( fdst0, v3constant, v3src0, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->dot( CVec3D * CVec3D[] )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->dot( fdst1, v3constant, v3src0, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( CMath::fabs( fdst0[i] - fdst1[i] ) > 1e-5f ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->dot( CVec3D * CVec3D[] ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );


	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->dot( fdst0, v3constant, v4src0, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->dot( CVec3D * CPlane[] )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->dot( fdst1, v3constant, v4src0, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( CMath::fabs( fdst0[i] - fdst1[i] ) > 1e-5f ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->dot( CVec3D * CPlane[] ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );


	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->dot( fdst0, v3constant, drawVerts, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->dot( CVec3D * CVertex[] )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->dot( fdst1, v3constant, drawVerts, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( CMath::fabs( fdst0[i] - fdst1[i] ) > 1e-5f ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->dot( CVec3D * CVertex[] ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );


	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->dot( fdst0, v4constant, v3src0, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->dot( CPlane * CVec3D[] )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->dot( fdst1, v4constant, v3src0, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( CMath::fabs( fdst0[i] - fdst1[i] ) > 1e-5f ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->dot( CPlane * CVec3D[] ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );


	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->dot( fdst0, v4constant, v4src0, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->dot( CPlane * CPlane[] )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->dot( fdst1, v4constant, v4src0, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( CMath::fabs( fdst0[i] - fdst1[i] ) > 1e-5f ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->dot( CPlane * CPlane[] ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );


	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->dot( fdst0, v4constant, drawVerts, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->dot( CPlane * CVertex[] )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->dot( fdst1, v4constant, drawVerts, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( CMath::fabs( fdst0[i] - fdst1[i] ) > 1e-5f ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->dot( CPlane * CVertex[] ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );


	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->dot( fdst0, v3src0, v3src1, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->dot( CVec3D[] * CVec3D[] )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->dot( fdst1, v3src0, v3src1, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( CMath::fabs( fdst0[i] - fdst1[i] ) > 1e-4f ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->dot( CVec3D[] * CVec3D[] ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );


	Debug::debug(Debug::math,__FUNCTION__) <<"====================================" <<endl;

	float dot1 = 0.0f, dot2 = 0.0f;
	for ( j = 0; j < 50 && j < COUNT; j++ ) {

		bestClocksGeneric = 0;
		for ( i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			p_generic->dot( dot1, fsrc0, fsrc1, j );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
		}
		CProcClock::printClocks( va( "generic->dot( float[%2d] * float[%2d] )", j, j ), 1, bestClocksGeneric );

		bestClocksSIMD = 0;
		for ( i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			p_simd->dot( dot2, fsrc0, fsrc1, j );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		}
		result = CMath::fabs( dot1 - dot2 ) < 1e-4f ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->dot( float[%2d] * float[%2d] ) %s", j, j, result ), 1, bestClocksSIMD, bestClocksGeneric );
	}
}

/*
============
TestCompare
============
*/
void TestCompare() {
	int i;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( float fsrc0[COUNT] );
	ALIGN16( sf_u8 bytedst[COUNT] );
	ALIGN16( sf_u8 bytedst2[COUNT] );
	const char *result;

	CRandom srnd( RANDOM_SEED );

	for ( i = 0; i < COUNT; i++ ) {
		fsrc0[i] = srnd.randomFloat() * 10.0f;
	}

	Debug::debug(Debug::math,__FUNCTION__) <<"====================================\n" <<endl;

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->cmpGT( bytedst, fsrc0, 0.0f, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->cmpGT( float[] >= float )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->cmpGT( bytedst2, fsrc0, 0.0f, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( bytedst[i] != bytedst2[i] ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->cmpGT( float[] >= float ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		memset( bytedst, 0, COUNT );
		StartRecordTimeLocal( start );
		p_generic->cmpGT( bytedst, 2, fsrc0, 0.0f, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->cmpGT( 2, float[] >= float )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		memset( bytedst2, 0, COUNT );
		StartRecordTimeLocal( start );
		p_simd->cmpGT( bytedst2, 2, fsrc0, 0.0f, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( bytedst[i] != bytedst2[i] ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->cmpGT( 2, float[] >= float ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );

	// ======================

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->cmpGE( bytedst, fsrc0, 0.0f, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->cmpGE( float[] >= float )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->cmpGE( bytedst2, fsrc0, 0.0f, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( bytedst[i] != bytedst2[i] ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->cmpGE( float[] >= float ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		memset( bytedst, 0, COUNT );
		StartRecordTimeLocal( start );
		p_generic->cmpGE( bytedst, 2, fsrc0, 0.0f, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->cmpGE( 2, float[] >= float )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		memset( bytedst2, 0, COUNT );
		StartRecordTimeLocal( start );
		p_simd->cmpGE( bytedst2, 2, fsrc0, 0.0f, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( bytedst[i] != bytedst2[i] ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->cmpGE( 2, float[] >= float ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );

	// ======================

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->cmpLT( bytedst, fsrc0, 0.0f, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->cmpLT( float[] >= float )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->cmpLT( bytedst2, fsrc0, 0.0f, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( bytedst[i] != bytedst2[i] ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->cmpLT( float[] >= float ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		memset( bytedst, 0, COUNT );
		StartRecordTimeLocal( start );
		p_generic->cmpLT( bytedst, 2, fsrc0, 0.0f, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->cmpLT( 2, float[] >= float )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		memset( bytedst2, 0, COUNT );
		StartRecordTimeLocal( start );
		p_simd->cmpLT( bytedst2, 2, fsrc0, 0.0f, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( bytedst[i] != bytedst2[i] ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->cmpLT( 2, float[] >= float ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );

	// ======================

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->cmpLE( bytedst, fsrc0, 0.0f, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->cmpLE( float[] >= float )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->cmpLE( bytedst2, fsrc0, 0.0f, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( bytedst[i] != bytedst2[i] ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->cmpLE( float[] >= float ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		memset( bytedst, 0, COUNT );
		StartRecordTimeLocal( start );
		p_generic->cmpLE( bytedst, 2, fsrc0, 0.0f, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->cmpLE( 2, float[] >= float )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		memset( bytedst2, 0, COUNT );
		StartRecordTimeLocal( start );
		p_simd->cmpLE( bytedst2, 2, fsrc0, 0.0f, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( bytedst[i] != bytedst2[i] ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->cmpLE( 2, float[] >= float ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );
}

/*
============
TestMinMax
============
*/
void TestMinMax() {
	int i;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( float fsrc0[COUNT] );
	ALIGN16( CVec2D v2src0[COUNT] );
	ALIGN16( CVec3D v3src0[COUNT] );
	ALIGN16( CVertex drawVerts[COUNT] );
	ALIGN16( int indexes[COUNT] );
	float min = 0.0f, max = 0.0f, min2 = 0.0f, max2 = 0.0f;
	CVec2D v2min, v2max, v2min2, v2max2;
	CVec3D vmin, vmax, vmin2, vmax2;
	const char *result;

	CRandom srnd( RANDOM_SEED );

	for ( i = 0; i < COUNT; i++ ) {
		fsrc0[i] = srnd.randomFloat() * 10.0f;
		v2src0[i][0] = srnd.randomFloat() * 10.0f;
		v2src0[i][1] = srnd.randomFloat() * 10.0f;
		v3src0[i][0] = srnd.randomFloat() * 10.0f;
		v3src0[i][1] = srnd.randomFloat() * 10.0f;
		v3src0[i][2] = srnd.randomFloat() * 10.0f;
		drawVerts[i].xyz = v3src0[i];
		indexes[i] = i;
	}

	Debug::debug(Debug::math,__FUNCTION__) <<"====================================\n"<<endl;

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		min = CMath::INFINITY_FLOAT;
		max = -CMath::INFINITY_FLOAT;
		StartRecordTimeLocal( start );
		p_generic->minMax( min, max, fsrc0, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->minMax( float[] )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->minMax( min2, max2, fsrc0, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	result = ( min == min2 && max == max2 ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->minMax( float[] ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->minMax( v2min, v2max, v2src0, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->minMax( CVec2D[] )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->minMax( v2min2, v2max2, v2src0, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	result = ( v2min == v2min2 && v2max == v2max2 ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->minMax( CVec2D[] ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->minMax( vmin, vmax, v3src0, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->minMax( CVec3D[] )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->minMax( vmin2, vmax2, v3src0, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	result = ( vmin == vmin2 && vmax == vmax2 ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->minMax( CVec3D[] ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->minMax( vmin, vmax, drawVerts, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->minMax( CVertex[] )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->minMax( vmin2, vmax2, drawVerts, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	result = ( vmin == vmin2 && vmax == vmax2 ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->minMax( CVertex[] ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->minMax( vmin, vmax, drawVerts, indexes, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->minMax( CVertex[], indexes[] )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->minMax( vmin2, vmax2, drawVerts, indexes, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	result = ( vmin == vmin2 && vmax == vmax2 ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->minMax( CVertex[], indexes[] ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );
}

/*
============
TestClamp
============
*/
void TestClamp() {
	int i;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( float fdst0[COUNT] );
	ALIGN16( float fdst1[COUNT] );
	ALIGN16( float fsrc0[COUNT] );
	const char *result;

	CRandom srnd( RANDOM_SEED );

	for ( i = 0; i < COUNT; i++ ) {
		fsrc0[i] = srnd.randomFloat() * 10.0f;
	}

	Debug::debug(Debug::math,__FUNCTION__) <<"====================================\n"<<endl;

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->clamp( fdst0, fsrc0, -1.0f, 1.0f, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->clamp( float[] )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->clamp( fdst1, fsrc0, -1.0f, 1.0f, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( fdst0[i] != fdst1[i] ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->clamp( float[] ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );


	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->clampMin( fdst0, fsrc0, -1.0f, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->clampMin( float[] )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->clampMin( fdst1, fsrc0, -1.0f, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( fdst0[i] != fdst1[i] ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->clampMin( float[] ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );


	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->clampMax( fdst0, fsrc0, 1.0f, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->clampMax( float[] )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->clampMax( fdst1, fsrc0, 1.0f, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( fdst0[i] != fdst1[i] ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->clampMax( float[] ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );
}

/*
============
TestMemcpy
============
*/
void TestMemcpy() {
	int i, j;
	sf_u8 test0[8192];
	sf_u8 test1[8192];

	CRandom random( RANDOM_SEED );

	Debug::debug(Debug::math,__FUNCTION__) <<"====================================\n" <<endl;

	for ( i = 5; i < 8192; i += 31 ) {
		for ( j = 0; j < i; j++ ) {
			test0[j] = random.randomInt( 255 );
		}
		p_simd->memCopy( test1, test0, 8192 );
		for ( j = 0; j < i; j++ ) {
			if ( test1[j] != test0[j] ) {
				Debug::debug(Debug::math,__FUNCTION__) << "   simd->memCopy() "<< S_COLOR_RED<<"X\n" <<endl;
				return;
			}
		}
	}
	Debug::debug(Debug::math,__FUNCTION__) << "   simd->memCopy() ok\n"<<endl;
}

/*
============
TestMemset
============
*/
void TestMemset() {
	int i, j, k;
	sf_u8 test[8192];

	for ( i = 0; i < 8192; i++ ) {
		test[i] = 0;
	}

	for ( i = 5; i < 8192; i += 31 ) {
		for ( j = -1; j <= 1; j++ ) {
			p_simd->memSet( test, j, i );
			for ( k = 0; k < i; k++ ) {
				if ( test[k] != (sf_u8)j ) {
					Debug::debug(Debug::math,__FUNCTION__) << "   simd->memSet() "<< S_COLOR_RED<< "X\n" <<endl;
					return;
				}
			}
		}
	}
	Debug::debug(Debug::math,__FUNCTION__) << "   simd->memSet() ok\n" <<endl;
}


/*
============
TestMatXMultiplyVecX
============
*/
void TestMatXMultiplyVecX() {
	int i, j;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	const char *result;
	CMatXD mat;
	CVecXD src(6);
	CVecXD dst(6);
	CVecXD tst(6);

	src[0] = 1.0f;
	src[1] = 2.0f;
	src[2] = 3.0f;
	src[3] = 4.0f;
	src[4] = 5.0f;
	src[5] = 6.0f;

	Debug::debug(Debug::math,__FUNCTION__) <<"================= NxN * Nx1 ===================\n"<<endl;

	for ( i = 1; i <= 6; i++ ) {
		mat.random( i, i, RANDOM_SEED, -10.0f, 10.0f );

		bestClocksGeneric = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			dst.toZero();
			StartRecordTimeLocal( start );
			p_generic->matX_MultiplyVecX( dst, mat, src );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
		}
		tst = dst;

		CProcClock::printClocks( va( "generic->matX_MultiplyVecX %dx%d*%dx1", i, i, i ), 1, bestClocksGeneric );

		bestClocksSIMD = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->matX_MultiplyVecX( dst, mat, src );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		}

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->matX_MultiplyVecX %dx%d*%dx1 %s", i, i, i, result ), 1, bestClocksSIMD, bestClocksGeneric );
	}

	Debug::debug(Debug::math,__FUNCTION__) <<"================= Nx6 * 6x1 ===================\n" <<endl;

	for ( i = 1; i <= 6; i++ ) {
		mat.random( i, 6, RANDOM_SEED, -10.0f, 10.0f );

		bestClocksGeneric = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			dst.toZero();
			StartRecordTimeLocal( start );
			p_generic->matX_MultiplyVecX( dst, mat, src );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
		}
		tst = dst;

		CProcClock::printClocks( va( "generic->matX_MultiplyVecX %dx6*6x1", i ), 1, bestClocksGeneric );

		bestClocksSIMD = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->matX_MultiplyVecX( dst, mat, src );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		}

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->matX_MultiplyVecX %dx6*6x1 %s", i, result ), 1, bestClocksSIMD, bestClocksGeneric );
	}

	Debug::debug(Debug::math,__FUNCTION__) <<"================= 6xN * Nx1 ===================\n"<<endl;

	for ( i = 1; i <= 6; i++ ) {
		mat.random( 6, i, RANDOM_SEED, -10.0f, 10.0f );

		bestClocksGeneric = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			dst.toZero();
			StartRecordTimeLocal( start );
			p_generic->matX_MultiplyVecX( dst, mat, src );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
		}
		tst = dst;

		CProcClock::printClocks( va( "generic->matX_MultiplyVecX 6x%d*%dx1", i, i ), 1, bestClocksGeneric );

		bestClocksSIMD = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			StartRecordTimeLocal( start );
			p_simd->matX_MultiplyVecX( dst, mat, src );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		}

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->matX_MultiplyVecX 6x%d*%dx1 %s", i, i, result ), 1, bestClocksSIMD, bestClocksGeneric );
	}
}

/*
============
TestMatXMultiplyAddVecX
============
*/
void TestMatXMultiplyAddVecX() {
	int i, j;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	const char *result;
	CMatXD mat;
	CVecXD src(6);
	CVecXD dst(6);
	CVecXD tst(6);

	src[0] = 1.0f;
	src[1] = 2.0f;
	src[2] = 3.0f;
	src[3] = 4.0f;
	src[4] = 5.0f;
	src[5] = 6.0f;

	Debug::debug(Debug::math,__FUNCTION__) <<"================= NxN * Nx1 ===================\n" <<endl;

	for ( i = 1; i <= 6; i++ ) {
		mat.random( i, i, RANDOM_SEED, -10.0f, 10.0f );

		bestClocksGeneric = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			dst.toZero();
			StartRecordTimeLocal( start );
			p_generic->matX_MultiplyAddVecX( dst, mat, src );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
		}
		tst = dst;

		CProcClock::printClocks( va( "generic->matX_MultiplyAddVecX %dx%d*%dx1", i, i, i ), 1, bestClocksGeneric );

		bestClocksSIMD = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->matX_MultiplyAddVecX( dst, mat, src );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		}

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->matX_MultiplyAddVecX %dx%d*%dx1 %s", i, i, i, result ), 1, bestClocksSIMD, bestClocksGeneric );
	}

	Debug::debug(Debug::math,__FUNCTION__) <<"================= Nx6 * 6x1 ===================\n" <<endl;

	for ( i = 1; i <= 6; i++ ) {
		mat.random( i, 6, RANDOM_SEED, -10.0f, 10.0f );

		bestClocksGeneric = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			dst.toZero();
			StartRecordTimeLocal( start );
			p_generic->matX_MultiplyAddVecX( dst, mat, src );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
		}
		tst = dst;

		CProcClock::printClocks( va( "generic->matX_MultiplyAddVecX %dx6*6x1", i ), 1, bestClocksGeneric );

		bestClocksSIMD = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->matX_MultiplyAddVecX( dst, mat, src );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		}

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->matX_MultiplyAddVecX %dx6*6x1 %s", i, result ), 1, bestClocksSIMD, bestClocksGeneric );
	}

	Debug::debug(Debug::math,__FUNCTION__) <<"================= 6xN * Nx1 ===================\n" <<endl;

	for ( i = 1; i <= 6; i++ ) {
		mat.random( 6, i, RANDOM_SEED, -10.0f, 10.0f );

		bestClocksGeneric = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			dst.toZero();
			StartRecordTimeLocal( start );
			p_generic->matX_MultiplyAddVecX( dst, mat, src );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
		}
		tst = dst;

		CProcClock::printClocks( va( "generic->matX_MultiplyAddVecX 6x%d*%dx1", i, i ), 1, bestClocksGeneric );

		bestClocksSIMD = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->matX_MultiplyAddVecX( dst, mat, src );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		}

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->matX_MultiplyAddVecX 6x%d*%dx1 %s", i, i, result ), 1, bestClocksSIMD, bestClocksGeneric );
	}
}

/*
============
TestMatXTransposeMultiplyVecX
============
*/
void TestMatXTransposeMultiplyVecX() {
	int i, j;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	const char *result;
	CMatXD mat;
	CVecXD src(6);
	CVecXD dst(6);
	CVecXD tst(6);

	src[0] = 1.0f;
	src[1] = 2.0f;
	src[2] = 3.0f;
	src[3] = 4.0f;
	src[4] = 5.0f;
	src[5] = 6.0f;

	Debug::debug(Debug::math,__FUNCTION__) <<"================= Nx6 * Nx1 ===================\n"<<endl;

	for ( i = 1; i <= 6; i++ ) {
		mat.random( i, 6, RANDOM_SEED, -10.0f, 10.0f );

		bestClocksGeneric = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			dst.toZero();
			StartRecordTimeLocal( start );
			p_generic->matX_TransposeMultiplyVecX( dst, mat, src );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
		}
		tst = dst;

		CProcClock::printClocks( va( "generic->matX_TransposeMulVecX %dx6*%dx1", i, i ), 1, bestClocksGeneric );

		bestClocksSIMD = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->matX_TransposeMultiplyVecX( dst, mat, src );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		}

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->matX_TransposeMulVecX %dx6*%dx1 %s", i, i, result ), 1, bestClocksSIMD, bestClocksGeneric );
	}

	Debug::debug(Debug::math,__FUNCTION__) <<"================= 6xN * 6x1 ===================\n" <<endl;

	for ( i = 1; i <= 6; i++ ) {
		mat.random( 6, i, RANDOM_SEED, -10.0f, 10.0f );

		bestClocksGeneric = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			dst.toZero();
			StartRecordTimeLocal( start );
			p_generic->matX_TransposeMultiplyVecX( dst, mat, src );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
		}
		tst = dst;

		CProcClock::printClocks( va( "generic->matX_TransposeMulVecX 6x%d*6x1", i ), 1, bestClocksGeneric );

		bestClocksSIMD = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->matX_TransposeMultiplyVecX( dst, mat, src );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		}

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->matX_TransposeMulVecX 6x%d*6x1 %s", i, result ), 1, bestClocksSIMD, bestClocksGeneric );
	}
}

/*
============
TestMatXTransposeMultiplyAddVecX
============
*/
void TestMatXTransposeMultiplyAddVecX() {
	int i, j;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	const char *result;
	CMatXD mat;
	CVecXD src(6);
	CVecXD dst(6);
	CVecXD tst(6);

	src[0] = 1.0f;
	src[1] = 2.0f;
	src[2] = 3.0f;
	src[3] = 4.0f;
	src[4] = 5.0f;
	src[5] = 6.0f;

	Debug::debug(Debug::math,__FUNCTION__) <<"================= Nx6 * Nx1 ===================\n" <<endl;

	for ( i = 1; i <= 6; i++ ) {
		mat.random( i, 6, RANDOM_SEED, -10.0f, 10.0f );

		bestClocksGeneric = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			dst.toZero();
			StartRecordTimeLocal( start );
			p_generic->matX_TransposeMultiplyAddVecX( dst, mat, src );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
		}
		tst = dst;

		CProcClock::printClocks( va( "generic->matX_TransposeMulAddVecX %dx6*%dx1", i, i ), 1, bestClocksGeneric );

		bestClocksSIMD = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->matX_TransposeMultiplyAddVecX( dst, mat, src );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		}

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->matX_TransposeMulAddVecX %dx6*%dx1 %s", i, i, result ), 1, bestClocksSIMD, bestClocksGeneric );
	}

	Debug::debug(Debug::math,__FUNCTION__) <<"================= 6xN * 6x1 ===================\n" <<endl;

	for ( i = 1; i <= 6; i++ ) {
		mat.random( 6, i, RANDOM_SEED, -10.0f, 10.0f );

		bestClocksGeneric = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			dst.toZero();
			StartRecordTimeLocal( start );
			p_generic->matX_TransposeMultiplyAddVecX( dst, mat, src );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
		}
		tst = dst;

		CProcClock::printClocks( va( "generic->matX_TransposeMulAddVecX 6x%d*6x1", i ), 1, bestClocksGeneric );

		bestClocksSIMD = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->matX_TransposeMultiplyAddVecX( dst, mat, src );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		}

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->matX_TransposeMulAddVecX 6x%d*6x1 %s", i, result ), 1, bestClocksSIMD, bestClocksGeneric );
	}
}

/*
============
TestMatXMultiplyMatX
============
*/
#define TEST_VALUE_RANGE			10.0f
#define	MATX_MATX_SIMD_EPSILON		1e-4f

void TestMatXMultiplyMatX() {
	int i, j;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	const char *result;
	CMatXD m1, m2, dst, tst;

	Debug::debug(Debug::math,__FUNCTION__) <<"================= NxN * Nx6 ===================" <<endl;

	// NxN * Nx6
	for ( i = 1; i <= 5; i++ ) {
		m1.random( i, i, RANDOM_SEED, -TEST_VALUE_RANGE, TEST_VALUE_RANGE );
		m2.random( i, 6, RANDOM_SEED, -TEST_VALUE_RANGE, TEST_VALUE_RANGE );
		dst.setSize( i, 6 );

		bestClocksGeneric = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			StartRecordTimeLocal( start );
			p_generic->matX_MultiplyMatX( dst, m1, m2 );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
		}
		tst = dst;

		CProcClock::printClocks( va( "generic->matX_MultiplyMatX %dx%d*%dx6", i, i, i ), 1, bestClocksGeneric );

		bestClocksSIMD = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			StartRecordTimeLocal( start );
			p_simd->matX_MultiplyMatX( dst, m1, m2 );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		}

		result = dst.compare( tst, MATX_MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->matX_MultiplyMatX %dx%d*%dx6 %s", i, i, i, result ), 1, bestClocksSIMD, bestClocksGeneric );
	}

	Debug::debug(Debug::math,__FUNCTION__) <<"================= 6xN * Nx6 ===================\n" <<endl;

	// 6xN * Nx6
	for ( i = 1; i <= 5; i++ ) {
		m1.random( 6, i, RANDOM_SEED, -TEST_VALUE_RANGE, TEST_VALUE_RANGE );
		m2.random( i, 6, RANDOM_SEED, -TEST_VALUE_RANGE, TEST_VALUE_RANGE );
		dst.setSize( 6, 6 );

		bestClocksGeneric = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			StartRecordTimeLocal( start );
			p_generic->matX_MultiplyMatX( dst, m1, m2 );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
		}
		tst = dst;

		CProcClock::printClocks( va( "generic->matX_MultiplyMatX 6x%d*%dx6", i, i ), 1, bestClocksGeneric );

		bestClocksSIMD = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			StartRecordTimeLocal( start );
			p_simd->matX_MultiplyMatX( dst, m1, m2 );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		}

		result = dst.compare( tst, MATX_MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->matX_MultiplyMatX 6x%d*%dx6 %s", i, i, result ), 1, bestClocksSIMD, bestClocksGeneric );
	}

	Debug::debug(Debug::math,__FUNCTION__) <<"================= Nx6 * 6xN ===================\n" <<endl;

	// Nx6 * 6xN
	for ( i = 1; i <= 5; i++ ) {
		m1.random( i, 6, RANDOM_SEED, -TEST_VALUE_RANGE, TEST_VALUE_RANGE );
		m2.random( 6, i, RANDOM_SEED, -TEST_VALUE_RANGE, TEST_VALUE_RANGE );
		dst.setSize( i, i );

		bestClocksGeneric = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			StartRecordTimeLocal( start );
			p_generic->matX_MultiplyMatX( dst, m1, m2 );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
		}
		tst = dst;

		CProcClock::printClocks( va( "generic->matX_MultiplyMatX %dx6*6x%d", i, i ), 1, bestClocksGeneric );

		bestClocksSIMD = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			StartRecordTimeLocal( start );
			p_simd->matX_MultiplyMatX( dst, m1, m2 );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		}

		result = dst.compare( tst, MATX_MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->matX_MultiplyMatX %dx6*6x%d %s", i, i, result ), 1, bestClocksSIMD, bestClocksGeneric );
	}

	Debug::debug(Debug::math,__FUNCTION__) <<"================= 6x6 * 6xN ===================\n" <<endl;

	// 6x6 * 6xN
	for ( i = 1; i <= 6; i++ ) {
		m1.random( 6, 6, RANDOM_SEED, -TEST_VALUE_RANGE, TEST_VALUE_RANGE );
		m2.random( 6, i, RANDOM_SEED, -TEST_VALUE_RANGE, TEST_VALUE_RANGE );
		dst.setSize( 6, i );

		bestClocksGeneric = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			StartRecordTimeLocal( start );
			p_generic->matX_MultiplyMatX( dst, m1, m2 );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
		}
		tst = dst;

		CProcClock::printClocks( va( "generic->matX_MultiplyMatX 6x6*6x%d", i ), 1, bestClocksGeneric );

		bestClocksSIMD = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			StartRecordTimeLocal( start );
			p_simd->matX_MultiplyMatX( dst, m1, m2 );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		}

		result = dst.compare( tst, MATX_MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->matX_MultiplyMatX 6x6*6x%d %s", i, result ), 1, bestClocksSIMD, bestClocksGeneric );
	}
}

/*
============
TestMatXTransposeMultiplyMatX
============
*/
void TestMatXTransposeMultiplyMatX() {
	int i, j;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	const char *result;
	CMatXD m1, m2, dst, tst;

	Debug::debug(Debug::math,__FUNCTION__) <<"================= Nx6 * NxN ===================\n" <<endl;

	// Nx6 * NxN
	for ( i = 1; i <= 5; i++ ) {
		m1.random( i, 6, RANDOM_SEED, -TEST_VALUE_RANGE, TEST_VALUE_RANGE );
		m2.random( i, i, RANDOM_SEED, -TEST_VALUE_RANGE, TEST_VALUE_RANGE );
		dst.setSize( 6, i );

		bestClocksGeneric = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			StartRecordTimeLocal( start );
			p_generic->matX_TransposeMultiplyMatX( dst, m1, m2 );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
		}
		tst = dst;

		CProcClock::printClocks( va( "generic->matX_TransMultiplyMatX %dx6*%dx%d", i, i, i ), 1, bestClocksGeneric );

		bestClocksSIMD = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			StartRecordTimeLocal( start );
			p_simd->matX_TransposeMultiplyMatX( dst, m1, m2 );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		}

		result = dst.compare( tst, MATX_MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->matX_TransMultiplyMatX %dx6*%dx%d %s", i, i, i, result ), 1, bestClocksSIMD, bestClocksGeneric );
	}

	Debug::debug(Debug::math,__FUNCTION__) <<"================= 6xN * 6x6 ===================\n" <<endl;

	// 6xN * 6x6
	for ( i = 1; i <= 6; i++ ) {
		m1.random( 6, i, RANDOM_SEED, -TEST_VALUE_RANGE, TEST_VALUE_RANGE );
		m2.random( 6, 6, RANDOM_SEED, -TEST_VALUE_RANGE, TEST_VALUE_RANGE );
		dst.setSize( i, 6 );

		bestClocksGeneric = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			StartRecordTimeLocal( start );
			p_generic->matX_TransposeMultiplyMatX( dst, m1, m2 );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
		}
		tst = dst;

		CProcClock::printClocks( va( "generic->matX_TransMultiplyMatX 6x%d*6x6", i ), 1, bestClocksGeneric );

		bestClocksSIMD = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			StartRecordTimeLocal( start );
			p_simd->matX_TransposeMultiplyMatX( dst, m1, m2 );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		}

		result = dst.compare( tst, MATX_MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->matX_TransMultiplyMatX 6x%d*6x6 %s", i, result ), 1, bestClocksSIMD, bestClocksGeneric );
	}
}

#define MATX_LTS_SIMD_EPSILON		1.0f
#define MATX_LTS_SOLVE_SIZE			100

/*
============
TestMatXLowerTriangularsolve
============
*/
void TestMatXLowerTriangularsolve() {
	int i, j;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	const char *result;
	CMatXD L;
	CVecXD x, b, tst;

	Debug::debug(Debug::math,__FUNCTION__) <<"====================================\n" <<endl;

	L.random( MATX_LTS_SOLVE_SIZE, MATX_LTS_SOLVE_SIZE, 0, -1.0f, 1.0f );
	x.setSize( MATX_LTS_SOLVE_SIZE );
	b.random( MATX_LTS_SOLVE_SIZE, 0, -1.0f, 1.0f );

	for ( i = 1; i < MATX_LTS_SOLVE_SIZE; i++ ) {

		x.zero( i );

		bestClocksGeneric = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			StartRecordTimeLocal( start );
			p_generic->matX_LowerTriangularsolve( L, x.toFloatPtr(), b.toFloatPtr(), i );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
		}
		tst = x;
		x.toZero();

		CProcClock::printClocks( va( "generic->matX_LowerTriangularsolve %dx%d", i, i ), 1, bestClocksGeneric );

		bestClocksSIMD = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			StartRecordTimeLocal( start );
			p_simd->matX_LowerTriangularsolve( L, x.toFloatPtr(), b.toFloatPtr(), i );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		}

		result = x.compare( tst, MATX_LTS_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->matX_LowerTriangularsolve %dx%d %s", i, i, result ), 1, bestClocksSIMD, bestClocksGeneric );
	}
}

/*
============
TestMatXLowerTriangularsolveTranspose
============
*/
void TestMatXLowerTriangularsolveTranspose() {
	int i, j;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	const char *result;
	CMatXD L;
	CVecXD x, b, tst;

	Debug::debug(Debug::math,__FUNCTION__) <<"====================================\n" <<endl;

	L.random( MATX_LTS_SOLVE_SIZE, MATX_LTS_SOLVE_SIZE, 0, -1.0f, 1.0f );
	x.setSize( MATX_LTS_SOLVE_SIZE );
	b.random( MATX_LTS_SOLVE_SIZE, 0, -1.0f, 1.0f );

	for ( i = 1; i < MATX_LTS_SOLVE_SIZE; i++ ) {

		x.zero( i );

		bestClocksGeneric = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			StartRecordTimeLocal( start );
			p_generic->matX_LowerTriangularsolveTranspose( L, x.toFloatPtr(), b.toFloatPtr(), i );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
		}
		tst = x;
		x.toZero();

		CProcClock::printClocks( va( "generic->matX_LowerTriangularsolveT %dx%d", i, i ), 1, bestClocksGeneric );

		bestClocksSIMD = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			StartRecordTimeLocal( start );
			p_simd->matX_LowerTriangularsolveTranspose( L, x.toFloatPtr(), b.toFloatPtr(), i );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		}

		result = x.compare( tst, MATX_LTS_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->matX_LowerTriangularsolveT %dx%d %s", i, i, result ), 1, bestClocksSIMD, bestClocksGeneric );
	}
}

#define MATX_ldlt_SIMD_EPSILON			0.1f
#define MATX_ldlt_FACTOR_SOLVE_SIZE		64

/*
============
TestMatXLDLTFactor
============
*/
void TestMatXLDLTFactor() {
	int i, j;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	const char *result;
	CMatXD src, original, mat1, mat2;
	CVecXD invDiag1, invDiag2;

	Debug::debug(Debug::math,__FUNCTION__) <<"====================================\n" <<endl;

	original.setSize( MATX_ldlt_FACTOR_SOLVE_SIZE, MATX_ldlt_FACTOR_SOLVE_SIZE );
	src.random( MATX_ldlt_FACTOR_SOLVE_SIZE, MATX_ldlt_FACTOR_SOLVE_SIZE, 0, -1.0f, 1.0f );
	src.transposeMultiply( original, src );

	for ( i = 1; i < MATX_ldlt_FACTOR_SOLVE_SIZE; i++ ) {

		bestClocksGeneric = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			mat1 = original;
			invDiag1.zero( MATX_ldlt_FACTOR_SOLVE_SIZE );
			StartRecordTimeLocal( start );
			p_generic->matX_LDLTFactor( mat1, invDiag1, i );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
		}

		CProcClock::printClocks( va( "generic->matX_LDLTFactor %dx%d", i, i ), 1, bestClocksGeneric );

		bestClocksSIMD = 0;
		for ( j = 0; j < NUMTESTS; j++ ) {
			mat2 = original;
			invDiag2.zero( MATX_ldlt_FACTOR_SOLVE_SIZE );
			StartRecordTimeLocal( start );
			p_simd->matX_LDLTFactor( mat2, invDiag2, i );
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		}

		result = mat1.compare( mat2, MATX_ldlt_SIMD_EPSILON ) && invDiag1.compare( invDiag2, MATX_ldlt_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->matX_LDLTFactor %dx%d %s", i, i, result ), 1, bestClocksSIMD, bestClocksGeneric );
	}
}

/*
============
TestblendJoints
============
*/
void TestblendJoints() {
	int i, j;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( CJointQuaternion baseJoints[COUNT] );
	ALIGN16( CJointQuaternion joints1[COUNT] );
	ALIGN16( CJointQuaternion joints2[COUNT] );
	ALIGN16( CJointQuaternion blendJoints[COUNT] );
	ALIGN16( int index[COUNT] );
	float lerp = 0.3f;
	const char *result;

	CRandom srnd( RANDOM_SEED );

	for ( i = 0; i < COUNT; i++ ) {
		CEulerAngles angles;
		angles[0] = srnd.randomFloat() * 180.0f;
		angles[1] = srnd.randomFloat() * 180.0f;
		angles[2] = srnd.randomFloat() * 180.0f;
		baseJoints[i].q = angles.toQuat();
		baseJoints[i].t[0] = srnd.randomFloat() * 10.0f;
		baseJoints[i].t[1] = srnd.randomFloat() * 10.0f;
		baseJoints[i].t[2] = srnd.randomFloat() * 10.0f;
		angles[0] = srnd.randomFloat() * 180.0f;
		angles[1] = srnd.randomFloat() * 180.0f;
		angles[2] = srnd.randomFloat() * 180.0f;
		blendJoints[i].q = angles.toQuat();
		blendJoints[i].t[0] = srnd.randomFloat() * 10.0f;
		blendJoints[i].t[1] = srnd.randomFloat() * 10.0f;
		blendJoints[i].t[2] = srnd.randomFloat() * 10.0f;
		index[i] = i;
	}

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		for ( j = 0; j < COUNT; j++ ) {
			joints1[j] = baseJoints[j];
		}
		StartRecordTimeLocal( start );
		p_generic->blendJoints( joints1, blendJoints, lerp, index, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->blendJoints()", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		for ( j = 0; j < COUNT; j++ ) {
			joints2[j] = baseJoints[j];
		}
		StartRecordTimeLocal( start );
		p_simd->blendJoints( joints2, blendJoints, lerp, index, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( !joints1[i].t.compare( joints2[i].t, 1e-3f ) ) {
			break;
		}
		if ( !joints1[i].q.compare( joints2[i].q, 1e-2f ) ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->blendJoints() %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );
}

/*
============
TestconvertJointQuatsToJointMats
============
*/
void TestconvertJointQuatsToJointMats() {
	int i;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( CJointQuaternion baseJoints[COUNT] );
	ALIGN16( CMatJoint3x4 joints1[COUNT] );
	ALIGN16( CMatJoint3x4 joints2[COUNT] );
	const char *result;

	CRandom srnd( RANDOM_SEED );

	for ( i = 0; i < COUNT; i++ ) {
		CEulerAngles angles;
		angles[0] = srnd.randomFloat() * 180.0f;
		angles[1] = srnd.randomFloat() * 180.0f;
		angles[2] = srnd.randomFloat() * 180.0f;
		baseJoints[i].q = angles.toQuat();
		baseJoints[i].t[0] = srnd.randomFloat() * 10.0f;
		baseJoints[i].t[1] = srnd.randomFloat() * 10.0f;
		baseJoints[i].t[2] = srnd.randomFloat() * 10.0f;
	}

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->convertJointQuatsToJointMats( joints1, baseJoints, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->convertJointQuatsToJointMats()", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->convertJointQuatsToJointMats( joints2, baseJoints, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( !joints1[i].compare( joints2[i], 1e-4f ) ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->convertJointQuatsToJointMats() %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );
}

/*
============
TestconvertJointMatsToJointQuats
============
*/
void TestconvertJointMatsToJointQuats() {
	int i;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( CMatJoint3x4 baseJoints[COUNT] );
	ALIGN16( CJointQuaternion joints1[COUNT] );
	ALIGN16( CJointQuaternion joints2[COUNT] );
	const char *result;

	CRandom srnd( RANDOM_SEED );

	for ( i = 0; i < COUNT; i++ ) {
		CEulerAngles angles;
		angles[0] = srnd.randomFloat() * 180.0f;
		angles[1] = srnd.randomFloat() * 180.0f;
		angles[2] = srnd.randomFloat() * 180.0f;
		baseJoints[i].setRotation( angles.toMat3() );
		CVec3D v;
		v[0] = srnd.randomFloat() * 10.0f;
		v[1] = srnd.randomFloat() * 10.0f;
		v[2] = srnd.randomFloat() * 10.0f;
		baseJoints[i].setTranslation( v );
	}

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->convertJointMatsToJointQuats( joints1, baseJoints, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->convertJointMatsToJointQuats()", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->convertJointMatsToJointQuats( joints2, baseJoints, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( !joints1[i].q.compare( joints2[i].q, 1e-4f ) ) {
			Debug::debug(Debug::math,__FUNCTION__) <<"convertJointMatsToJointQuats: broken q: "<< i <<endl;
			break;
		}
		if ( !joints1[i].t.compare( joints2[i].t, 1e-4f ) ) {
			Debug::debug(Debug::math,__FUNCTION__) <<"convertJointMatsToJointQuats: broken t: " << i <<endl;
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->convertJointMatsToJointQuats() %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );
}

/*
============
TesttransformJoints
============
*/
void TesttransformJoints() {
	int i, j;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( CMatJoint3x4 joints[COUNT+1] );
	ALIGN16( CMatJoint3x4 joints1[COUNT+1] );
	ALIGN16( CMatJoint3x4 joints2[COUNT+1] );
	ALIGN16( int parents[COUNT+1] );
	const char *result;

	CRandom srnd( RANDOM_SEED );

	for ( i = 0; i <= COUNT; i++ ) {
		CEulerAngles angles;
		angles[0] = srnd.randomFloat() * 180.0f;
		angles[1] = srnd.randomFloat() * 180.0f;
		angles[2] = srnd.randomFloat() * 180.0f;
		joints[i].setRotation( angles.toMat3() );
		CVec3D v;
		v[0] = srnd.randomFloat() * 2.0f;
		v[1] = srnd.randomFloat() * 2.0f;
		v[2] = srnd.randomFloat() * 2.0f;
		joints[i].setTranslation( v );
		parents[i] = i - 1;
	}

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		for ( j = 0; j <= COUNT; j++ ) {
			joints1[j] = joints[j];
		}
		StartRecordTimeLocal( start );
		p_generic->transformJoints( joints1, parents, 1, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->transformJoints()", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		for ( j = 0; j <= COUNT; j++ ) {
			joints2[j] = joints[j];
		}
		StartRecordTimeLocal( start );
		p_simd->transformJoints( joints2, parents, 1, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( !joints1[i+1].compare( joints2[i+1], 1e-4f ) ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->transformJoints() %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );
}

/*
============
TestuntransformJoints
============
*/
void TestuntransformJoints() {
	int i, j;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( CMatJoint3x4 joints[COUNT+1] );
	ALIGN16( CMatJoint3x4 joints1[COUNT+1] );
	ALIGN16( CMatJoint3x4 joints2[COUNT+1] );
	ALIGN16( int parents[COUNT+1] );
	const char *result;

	CRandom srnd( RANDOM_SEED );

	for ( i = 0; i <= COUNT; i++ ) {
		CEulerAngles angles;
		angles[0] = srnd.randomFloat() * 180.0f;
		angles[1] = srnd.randomFloat() * 180.0f;
		angles[2] = srnd.randomFloat() * 180.0f;
		joints[i].setRotation( angles.toMat3() );
		CVec3D v;
		v[0] = srnd.randomFloat() * 2.0f;
		v[1] = srnd.randomFloat() * 2.0f;
		v[2] = srnd.randomFloat() * 2.0f;
		joints[i].setTranslation( v );
		parents[i] = i - 1;
	}

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		for ( j = 0; j <= COUNT; j++ ) {
			joints1[j] = joints[j];
		}
		StartRecordTimeLocal( start );
		p_generic->untransformJoints( joints1, parents, 1, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->untransformJoints()", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		for ( j = 0; j <= COUNT; j++ ) {
			joints2[j] = joints[j];
		}
		StartRecordTimeLocal( start );
		p_simd->untransformJoints( joints2, parents, 1, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( !joints1[i+1].compare( joints2[i+1], 1e-4f ) ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->untransformJoints() %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );
}

/*
============
TesttransformVerts
============
*/
#define NUMJOINTS	64
#define NUMVERTS	COUNT/2
void TesttransformVerts() {
	int i;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( CVertex drawVerts1[NUMVERTS] );
	ALIGN16( CVertex drawVerts2[NUMVERTS] );
	ALIGN16( CMatJoint3x4 joints[NUMJOINTS] );
	ALIGN16( CVec4D weights[COUNT] );
	ALIGN16( int weightIndex[COUNT*2] );
	const char *result;

	CRandom srnd( RANDOM_SEED );

	for ( i = 0; i < NUMJOINTS; i++ ) {
		CEulerAngles angles;
		angles[0] = srnd.randomFloat() * 180.0f;
		angles[1] = srnd.randomFloat() * 180.0f;
		angles[2] = srnd.randomFloat() * 180.0f;
		joints[i].setRotation( angles.toMat3() );
		CVec3D v;
		v[0] = srnd.randomFloat() * 2.0f;
		v[1] = srnd.randomFloat() * 2.0f;
		v[2] = srnd.randomFloat() * 2.0f;
		joints[i].setTranslation( v );
	}

	for ( i = 0; i < COUNT; i++ ) {
		weights[i][0] = srnd.randomFloat() * 2.0f;
		weights[i][1] = srnd.randomFloat() * 2.0f;
		weights[i][2] = srnd.randomFloat() * 2.0f;
		weights[i][3] = srnd.randomFloat();
		weightIndex[i*2+0] = ( i * NUMJOINTS / COUNT ) * sizeof( CMatJoint3x4 );
		weightIndex[i*2+1] = i & 1;
	}

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->transformVerts( drawVerts1, NUMVERTS, joints, weights, weightIndex, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->transformVerts()", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->transformVerts( drawVerts2, NUMVERTS, joints, weights, weightIndex, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < NUMVERTS; i++ ) {
		if ( !drawVerts1[i].xyz.compare( drawVerts2[i].xyz, 0.5f ) ) {
			break;
		}
	}
	result = ( i >= NUMVERTS ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->transformVerts() %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );
}

/*
============
TesttracePointCull
============
*/
void TesttracePointCull() {
	int i, j;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( CPlane planes[4] );
	ALIGN16( CVertex drawVerts[COUNT] );
	ALIGN16( sf_u8 cullBits1[COUNT] );
	ALIGN16( sf_u8 cullBits2[COUNT] );
	sf_u8 totalOr1 = 0, totalOr2 = 0;
	const char *result;

	CRandom srnd( RANDOM_SEED );

	planes[0].setNormal( CVec3D(  1,  0, 0 ) );
	planes[1].setNormal( CVec3D( -1,  0, 0 ) );
	planes[2].setNormal( CVec3D(  0,  1, 0 ) );
	planes[3].setNormal( CVec3D(  0, -1, 0 ) );
	planes[0][3] = -5.3f;
	planes[1][3] = 5.3f;
	planes[2][3] = -3.4f;
	planes[3][3] = 3.4f;

	for ( i = 0; i < COUNT; i++ ) {
		for ( j = 0; j < 3; j++ ) {
			drawVerts[i].xyz[j] = srnd.randomFloat() * 10.0f;
		}
	}

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->tracePointCull( cullBits1, totalOr1, 0.0f, planes, drawVerts, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->tracePointCull()", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->tracePointCull( cullBits2, totalOr2, 0.0f, planes, drawVerts, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( cullBits1[i] != cullBits2[i] ) {
			break;
		}
	}
	result = ( i >= COUNT && totalOr1 == totalOr2 ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->tracePointCull() %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );
}

/*
============
TestdecalPointCull
============
*/
void TestdecalPointCull() {
	int i, j;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( CPlane planes[6] );
	ALIGN16( CVertex drawVerts[COUNT] );
	ALIGN16( sf_u8 cullBits1[COUNT] );
	ALIGN16( sf_u8 cullBits2[COUNT] );
	const char *result;

	CRandom srnd( RANDOM_SEED );

	planes[0].setNormal( CVec3D(  1,  0,  0 ) );
	planes[1].setNormal( CVec3D( -1,  0,  0 ) );
	planes[2].setNormal( CVec3D(  0,  1,  0 ) );
	planes[3].setNormal( CVec3D(  0, -1,  0 ) );
	planes[4].setNormal( CVec3D(  0,  0,  1 ) );
	planes[5].setNormal( CVec3D(  0,  0, -1 ) );
	planes[0][3] = -5.3f;
	planes[1][3] = 5.3f;
	planes[2][3] = -4.4f;
	planes[3][3] = 4.4f;
	planes[4][3] = -3.5f;
	planes[5][3] = 3.5f;

	for ( i = 0; i < COUNT; i++ ) {
		for ( j = 0; j < 3; j++ ) {
			drawVerts[i].xyz[j] = srnd.randomFloat() * 10.0f;
		}
	}

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->decalPointCull( cullBits1, planes, drawVerts, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->decalPointCull()", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->decalPointCull( cullBits2, planes, drawVerts, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( cullBits1[i] != cullBits2[i] ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->decalPointCull() %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );
}

/*
============
TestoverlayPointCull
============
*/
void TestoverlayPointCull() {
	int i, j;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( CPlane planes[2] );
	ALIGN16( CVertex drawVerts[COUNT] );
	ALIGN16( sf_u8 cullBits1[COUNT] );
	ALIGN16( sf_u8 cullBits2[COUNT] );
	ALIGN16( CVec2D texCoords1[COUNT] );
	ALIGN16( CVec2D texCoords2[COUNT] );
	const char *result;

	CRandom srnd( RANDOM_SEED );

	planes[0].setNormal( CVec3D( 0.3f, 0.2f, 0.9f ) );
	planes[1].setNormal( CVec3D( 0.9f, 0.2f, 0.3f ) );
	planes[0][3] = -5.3f;
	planes[1][3] = -4.3f;

	for ( i = 0; i < COUNT; i++ ) {
		for ( j = 0; j < 3; j++ ) {
			drawVerts[i].xyz[j] = srnd.randomFloat() * 10.0f;
		}
	}

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->overlayPointCull( cullBits1, texCoords1, planes, drawVerts, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->overlayPointCull()", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->overlayPointCull( cullBits2, texCoords2, planes, drawVerts, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( cullBits1[i] != cullBits2[i] ) {
			break;
		}
		if ( !texCoords1[i].compare( texCoords2[i], 1e-4f ) ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->overlayPointCull() %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );
}

/*
============
TestderiveTriPlanes
============
*/
void TestderiveTriPlanes() {
	int i, j;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( CVertex drawVerts1[COUNT] );
	ALIGN16( CVertex drawVerts2[COUNT] );
	ALIGN16( GEO::CPlane planes1[COUNT] );
	ALIGN16( GEO::CPlane planes2[COUNT] );
	ALIGN16( int indexes[COUNT*3] );
	const char *result;

	CRandom srnd( RANDOM_SEED );

	for ( i = 0; i < COUNT; i++ ) {
		for ( j = 0; j < 3; j++ ) {
			drawVerts1[i].xyz[j] = srnd.randomFloat() * 10.0f;
		}
		for ( j = 0; j < 2; j++ ) {
			drawVerts1[i].st[j] = srnd.randomFloat();
		}
		drawVerts2[i] = drawVerts1[i];
	}

	for ( i = 0; i < COUNT; i++ ) {
		indexes[i*3+0] = ( i + 0 ) % COUNT;
		indexes[i*3+1] = ( i + 1 ) % COUNT;
		indexes[i*3+2] = ( i + 2 ) % COUNT;
	}

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->deriveTriPlanes( planes1, drawVerts1, COUNT, indexes, COUNT*3 );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->deriveTriPlanes()", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->deriveTriPlanes( planes2, drawVerts2, COUNT, indexes, COUNT*3 );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( !planes1[i].compare( planes2[i], 1e-1f, 1e-1f ) ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->deriveTriPlanes() %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );
}

/*
============
Test_deriveTangents
============
*/
void Test_deriveTangents() {
	int i, j;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( CVertex drawVerts1[COUNT] );
	ALIGN16( CVertex drawVerts2[COUNT] );
	ALIGN16( SMF::GEO::CPlane planes1[COUNT] );
	ALIGN16( SMF::GEO::CPlane planes2[COUNT] );
	ALIGN16( int indexes[COUNT*3] );
	const char *result;

	CRandom srnd( RANDOM_SEED );

	for ( i = 0; i < COUNT; i++ ) {
		for ( j = 0; j < 3; j++ ) {
			drawVerts1[i].xyz[j] = srnd.randomFloat() * 10.0f;
		}
		for ( j = 0; j < 2; j++ ) {
			drawVerts1[i].st[j] = srnd.randomFloat();
		}
		drawVerts2[i] = drawVerts1[i];
	}

	for ( i = 0; i < COUNT; i++ ) {
		indexes[i*3+0] = ( i + 0 ) % COUNT;
		indexes[i*3+1] = ( i + 1 ) % COUNT;
		indexes[i*3+2] = ( i + 2 ) % COUNT;
	}

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic-> deriveTangents( planes1, drawVerts1, COUNT, indexes, COUNT*3 );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic-> deriveTangents()", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd-> deriveTangents( planes2, drawVerts2, COUNT, indexes, COUNT*3 );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		CVec3D v1, v2;

		v1 = drawVerts1[i].normal;
		v1.toNormal();
		v2 = drawVerts2[i].normal;
		v2.toNormal();
		if ( !v1.compare( v2, 1e-1f ) ) {
			Debug::debug(Debug::math,__FUNCTION__) <<" deriveTangents: broken at tangent0: "<< i<< " -- expecting "<< v1.toString()<<" got: " <<v2.toString()<<endl;
			break;
		}
		v1 = drawVerts1[i].tangents[0];
		v1.toNormal();
		v2 = drawVerts2[i].tangents[0];
		v2.toNormal();
		if ( !v1.compare( v2, 1e-1f ) ) {
			Debug::debug(Debug::math,__FUNCTION__) <<" deriveTangents: broken at tangent0: "<< i<< " -- expecting "<< v1.toString()<<" got: " <<v2.toString()<<endl;
			break;
		}
		v1 = drawVerts1[i].tangents[1];
		v1.toNormal();
		v2 = drawVerts2[i].tangents[1];
		v2.toNormal();
		if ( !v1.compare( v2, 1e-1f ) ) {
			Debug::debug(Debug::math,__FUNCTION__) <<" deriveTangents: broken at tangent1: "<< i<<" -- expecting " << v1.toString()<<" got: "<< v2.toString()<<endl;
			break;
		}
		if ( !planes1[i].compare( planes2[i], 1e-1f, 1e-1f ) ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd-> deriveTangents() %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );
}

/*
============
TestDeriveUnsmoothedTangents
============
*/
/*
void TestDeriveUnsmoothedTangents() {
	int i, j;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( CVertex drawVerts1[COUNT] );
	ALIGN16( CVertex drawVerts2[COUNT] );
	ALIGN16( dominantTri_s dominantTris[COUNT] );
	const char *result;

	CRandom srnd( RANDOM_SEED );

	for ( i = 0; i < COUNT; i++ ) {
		for ( j = 0; j < 3; j++ ) {
			drawVerts1[i].xyz[j] = srnd.randomFloat() * 10.0f;
		}
		for ( j = 0; j < 2; j++ ) {
			drawVerts1[i].st[j] = srnd.randomFloat();
		}
		drawVerts2[i] = drawVerts1[i];

		dominantTris[i].v2 = ( i + 1 + srnd.randomInt( 8 ) ) % COUNT;
		dominantTris[i].v3 = ( i + 9 + srnd.randomInt( 8 ) ) % COUNT;
		dominantTris[i].normalizationScale[0] = srnd.randomFloat();
		dominantTris[i].normalizationScale[1] = srnd.randomFloat();
		dominantTris[i].normalizationScale[2] = srnd.randomFloat();
	}

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->deriveUnsmoothedTangents( drawVerts1, dominantTris, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->deriveUnsmoothedTangents()", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->deriveUnsmoothedTangents( drawVerts2, dominantTris, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		CVec3D v1, v2;

		v1 = drawVerts1[i].normal;
		v1.toNormal();
		v2 = drawVerts2[i].normal;
		v2.toNormal();
		if ( !v1.compare( v2, 1e-1f ) ) {
			break;
		}
		v1 = drawVerts1[i].tangents[0];
		v1.toNormal();
		v2 = drawVerts2[i].tangents[0];
		v2.toNormal();
		if ( !v1.compare( v2, 1e-1f ) ) {
			break;
		}
		v1 = drawVerts1[i].tangents[1];
		v1.toNormal();
		v2 = drawVerts2[i].tangents[1];
		v2.toNormal();
		if ( !v1.compare( v2, 1e-1f ) ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->deriveUnsmoothedTangents() %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );
}
*/
/*
============
Test_normalizeTangents()
============
*/
void Test_normalizeTangents() {
	int i, j;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( CVertex drawVerts1[COUNT] );
	ALIGN16( CVertex drawVerts2[COUNT] );
	const char *result;

	CRandom srnd( RANDOM_SEED );

	for ( i = 0; i < COUNT; i++ ) {
		for ( j = 0; j < 3; j++ ) {
			drawVerts1[i].normal[j] = srnd.randomFloat() * 10.0f;
			drawVerts1[i].tangents[0][j] = srnd.randomFloat() * 10.0f;
			drawVerts1[i].tangents[1][j] = srnd.randomFloat() * 10.0f;
		}
		drawVerts2[i] = drawVerts1[i];
	}

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic-> normalizeTangents( drawVerts1, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic-> normalizeTangents()", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd-> normalizeTangents( drawVerts2, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( !drawVerts1[i].normal.compare( drawVerts2[i].normal, 1e-2f ) ) {
			break;
		}
		if ( !drawVerts1[i].tangents[0].compare( drawVerts2[i].tangents[0], 1e-2f ) ) {
			break;
		}
		if ( !drawVerts1[i].tangents[1].compare( drawVerts2[i].tangents[1], 1e-2f ) ) {
			break;
		}

		// since we're doing a lot of unaligned work, added this check to
		// make sure xyz wasn't getting overwritten
		if ( !drawVerts1[i].xyz.compare( drawVerts2[i].xyz, 1e-2f ) ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd-> normalizeTangents() %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );
}

/*
============
TestGetTextureSpaceLightVectors
============
*/
void TestGetTextureSpaceLightVectors() {
	int i, j;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( CVertex drawVerts[COUNT] );
	ALIGN16( CVec4D texCoords1[COUNT] );
	ALIGN16( CVec4D texCoords2[COUNT] );
	ALIGN16( int indexes[COUNT*3] );
	ALIGN16( CVec3D lightVectors1[COUNT] );
	ALIGN16( CVec3D lightVectors2[COUNT] );
	CVec3D lightOrigin;
	const char *result;

	CRandom srnd( RANDOM_SEED );

	for ( i = 0; i < COUNT; i++ ) {
		for ( j = 0; j < 3; j++ ) {
			drawVerts[i].xyz[j] = srnd.randomFloat() * 100.0f;
			drawVerts[i].normal[j] = srnd.randomFloat();
			drawVerts[i].tangents[0][j] = srnd.randomFloat();
			drawVerts[i].tangents[1][j] = srnd.randomFloat();
		}
	}

	for ( i = 0; i < COUNT; i++ ) {
		indexes[i*3+0] = ( i + 0 ) % COUNT;
		indexes[i*3+1] = ( i + 1 ) % COUNT;
		indexes[i*3+2] = ( i + 2 ) % COUNT;
	}

	lightOrigin[0] = srnd.randomFloat() * 100.0f;
	lightOrigin[1] = srnd.randomFloat() * 100.0f;
	lightOrigin[2] = srnd.randomFloat() * 100.0f;

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic-> createTextureSpaceLightVectors( lightVectors1, lightOrigin, drawVerts, COUNT, indexes, COUNT*3 );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic-> createTextureSpaceLightVectors()", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd-> createTextureSpaceLightVectors( lightVectors2, lightOrigin, drawVerts, COUNT, indexes, COUNT*3 );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( !lightVectors1[i].compare( lightVectors2[i], 1e-4f ) ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd-> createTextureSpaceLightVectors() %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );
}

/*
============
TestGetSpecularTextureCoords
============
*/
void TestGetSpecularTextureCoords() {
	int i, j;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( CVertex drawVerts[COUNT] );
	ALIGN16( CVec4D texCoords1[COUNT] );
	ALIGN16( CVec4D texCoords2[COUNT] );
	ALIGN16( int indexes[COUNT*3] );
	ALIGN16( CVec3D lightVectors1[COUNT] );
	ALIGN16( CVec3D lightVectors2[COUNT] );
	CVec3D lightOrigin, viewOrigin;
	const char *result;

	CRandom srnd( RANDOM_SEED );

	for ( i = 0; i < COUNT; i++ ) {
		for ( j = 0; j < 3; j++ ) {
			drawVerts[i].xyz[j] = srnd.randomFloat() * 100.0f;
			drawVerts[i].normal[j] = srnd.randomFloat();
			drawVerts[i].tangents[0][j] = srnd.randomFloat();
			drawVerts[i].tangents[1][j] = srnd.randomFloat();
		}
	}

	for ( i = 0; i < COUNT; i++ ) {
		indexes[i*3+0] = ( i + 0 ) % COUNT;
		indexes[i*3+1] = ( i + 1 ) % COUNT;
		indexes[i*3+2] = ( i + 2 ) % COUNT;
	}

	lightOrigin[0] = srnd.randomFloat() * 100.0f;
	lightOrigin[1] = srnd.randomFloat() * 100.0f;
	lightOrigin[2] = srnd.randomFloat() * 100.0f;
	viewOrigin[0] = srnd.randomFloat() * 100.0f;
	viewOrigin[1] = srnd.randomFloat() * 100.0f;
	viewOrigin[2] = srnd.randomFloat() * 100.0f;

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic-> createSpecularTextureCoords( texCoords1, lightOrigin, viewOrigin, drawVerts, COUNT, indexes, COUNT*3 );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic-> createSpecularTextureCoords()", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd-> createSpecularTextureCoords( texCoords2, lightOrigin, viewOrigin, drawVerts, COUNT, indexes, COUNT*3 );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( !texCoords1[i].compare( texCoords2[i], 1e-2f ) ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd-> createSpecularTextureCoords() %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );
}

/*
============
Test_createShadowCache
============
*/
void Test_createShadowCache() {
	int i, j;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( CVertex drawVerts[COUNT] );
	ALIGN16( CVec4D vertexCache1[COUNT*2] );
	ALIGN16( CVec4D vertexCache2[COUNT*2] );
	ALIGN16( int originalVertRemap[COUNT] );
	ALIGN16( int vertRemap1[COUNT] );
	ALIGN16( int vertRemap2[COUNT] );
	ALIGN16( CVec3D lightOrigin );
	int numVerts1 = 0, numVerts2 = 0;
	const char *result;

	CRandom srnd( RANDOM_SEED );

	for ( i = 0; i < COUNT; i++ ) {
		drawVerts[i].xyz[0] = srnd.randomFloat() * 100.0f;
		drawVerts[i].xyz[1] = srnd.randomFloat() * 100.0f;
		drawVerts[i].xyz[2] = srnd.randomFloat() * 100.0f;
		originalVertRemap[i] = ( srnd.randomFloat() > 0.0f ) ? -1 : 0;
	}
	lightOrigin[0] = srnd.randomFloat() * 100.0f;
	lightOrigin[1] = srnd.randomFloat() * 100.0f;
	lightOrigin[2] = srnd.randomFloat() * 100.0f;

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		for ( j = 0; j < COUNT; j++ ) {
			vertRemap1[j] = originalVertRemap[j];
		}
		StartRecordTimeLocal( start );
		numVerts1 =p_generic-> createShadowCache( vertexCache1, vertRemap1, lightOrigin, drawVerts, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic-> createShadowCache()", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		for ( j = 0; j < COUNT; j++ ) {
			vertRemap2[j] = originalVertRemap[j];
		}
		StartRecordTimeLocal( start );
		numVerts2 = p_simd-> createShadowCache( vertexCache2, vertRemap2, lightOrigin, drawVerts, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( i < ( numVerts1 / 2 ) ) {
			if ( !vertexCache1[i*2+0].compare( vertexCache2[i*2+0], 1e-2f ) ) {
				break;
			}
			if ( !vertexCache1[i*2+1].compare( vertexCache2[i*2+1], 1e-2f ) ) {
				break;
			}
		}
		if ( vertRemap1[i] != vertRemap2[i] ) {
			break;
		}
	}

	result = ( i >= COUNT && numVerts1 == numVerts2 ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd-> createShadowCache() %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_generic-> createVertexProgramShadowCache( vertexCache1, drawVerts, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic-> createVertexProgramShadowCache()", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		p_simd-> createVertexProgramShadowCache( vertexCache2, drawVerts, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( !vertexCache1[i*2+0].compare( vertexCache2[i*2+0], 1e-2f ) ) {
			break;
		}
		if ( !vertexCache1[i*2+1].compare( vertexCache2[i*2+1], 1e-2f ) ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd-> createVertexProgramShadowCache() %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );
}

/*
============
TestSoundUpSampling
============
*/
#define SOUND_UPSAMPLE_EPSILON		1.0f

void TestSoundUpSampling() {
	int i;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( short pcm[MIXBUFFER_SAMPLES*2] );
	ALIGN16( float ogg0[MIXBUFFER_SAMPLES*2] );
	ALIGN16( float ogg1[MIXBUFFER_SAMPLES*2] );
	ALIGN16( float samples1[MIXBUFFER_SAMPLES*2] );
	ALIGN16( float samples2[MIXBUFFER_SAMPLES*2] );
	float *ogg[2];
	int kHz, numSpeakers;
	const char *result;

	CRandom srnd( RANDOM_SEED );

	for ( i = 0; i < MIXBUFFER_SAMPLES*2; i++ ) {
		pcm[i] = srnd.randomInt( (1<<16) ) - (1<<15);
		ogg0[i] = srnd.randomFloat();
		ogg1[i] = srnd.randomFloat();
	}

	ogg[0] = ogg0;
	ogg[1] = ogg1;

	for ( numSpeakers = 1; numSpeakers <= 2; numSpeakers++ ) {

		for ( kHz = 11025; kHz <= 44100; kHz *= 2 ) {
			bestClocksGeneric = 0;
			for ( i = 0; i < NUMTESTS; i++ ) {
				StartRecordTimeLocal( start );
				p_generic-> upSamplePCMTo44kHz( samples1, pcm, MIXBUFFER_SAMPLES*numSpeakers*kHz/44100, kHz, numSpeakers );
				StopRecordTimeLocal( end );
				GetBest( start, end, bestClocksGeneric );
			}
			CProcClock::printClocks( va( "generic-> upSamplePCMTo44kHz( %d, %d )", kHz, numSpeakers ), MIXBUFFER_SAMPLES*numSpeakers*kHz/44100, bestClocksGeneric );

			bestClocksSIMD = 0;
			for ( i = 0; i < NUMTESTS; i++ ) {
				StartRecordTimeLocal( start );
				p_simd-> upSamplePCMTo44kHz( samples2, pcm, MIXBUFFER_SAMPLES*numSpeakers*kHz/44100, kHz, numSpeakers );
				StopRecordTimeLocal( end );
				GetBest( start, end, bestClocksSIMD );
			}

			for ( i = 0; i < MIXBUFFER_SAMPLES*numSpeakers; i++ ) {
				if ( CMath::fabs( samples1[i] - samples2[i] ) > SOUND_UPSAMPLE_EPSILON ) {
					break;
				}
			}
			result = ( i >= MIXBUFFER_SAMPLES*numSpeakers ) ? "ok" : S_COLOR_RED"X";
			CProcClock::printClocks( va( "   simd-> upSamplePCMTo44kHz( %d, %d ) %s", kHz, numSpeakers, result ), MIXBUFFER_SAMPLES*numSpeakers*kHz/44100, bestClocksSIMD, bestClocksGeneric );
		}
	}

	for ( numSpeakers = 1; numSpeakers <= 2; numSpeakers++ ) {

		for ( kHz = 11025; kHz <= 44100; kHz *= 2 ) {
			bestClocksGeneric = 0;
			for ( i = 0; i < NUMTESTS; i++ ) {
				StartRecordTimeLocal( start );
				p_generic-> upSampleOGGTo44kHz( samples1, ogg, MIXBUFFER_SAMPLES*numSpeakers*kHz/44100, kHz, numSpeakers );
				StopRecordTimeLocal( end );
				GetBest( start, end, bestClocksGeneric );
			}
			CProcClock::printClocks( va( "generic-> upSampleOGGTo44kHz( %d, %d )", kHz, numSpeakers ), MIXBUFFER_SAMPLES*numSpeakers*kHz/44100, bestClocksGeneric );

			bestClocksSIMD = 0;
			for ( i = 0; i < NUMTESTS; i++ ) {
				StartRecordTimeLocal( start );
				p_simd-> upSampleOGGTo44kHz( samples2, ogg, MIXBUFFER_SAMPLES*numSpeakers*kHz/44100, kHz, numSpeakers );
				StopRecordTimeLocal( end );
				GetBest( start, end, bestClocksSIMD );
			}

			for ( i = 0; i < MIXBUFFER_SAMPLES*numSpeakers; i++ ) {
				if ( CMath::fabs( samples1[i] - samples2[i] ) > SOUND_UPSAMPLE_EPSILON ) {
					break;
				}
			}
			result = ( i >= MIXBUFFER_SAMPLES ) ? "ok" : S_COLOR_RED"X";
			CProcClock::printClocks( va( "   simd-> upSampleOGGTo44kHz( %d, %d ) %s", kHz, numSpeakers, result ), MIXBUFFER_SAMPLES*numSpeakers*kHz/44100, bestClocksSIMD, bestClocksGeneric );
		}
	}
}

/*
============
TestSoundMixing
============
*/
#define SOUND_MIX_EPSILON		2.0f

void TestSoundMixing() {
	int i, j;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( float origMixBuffer[MIXBUFFER_SAMPLES*6] );
	ALIGN16( float mixBuffer1[MIXBUFFER_SAMPLES*6] );
	ALIGN16( float mixBuffer2[MIXBUFFER_SAMPLES*6] );
	ALIGN16( float samples[MIXBUFFER_SAMPLES*6] );
	ALIGN16( short outSamples1[MIXBUFFER_SAMPLES*6] );
	ALIGN16( short outSamples2[MIXBUFFER_SAMPLES*6] );
	float lastV[6];
	float currentV[6];
	const char *result;

	CRandom srnd( RANDOM_SEED );

	for ( i = 0; i < 6; i++ ) {
		lastV[i] = srnd.randomFloat();
		currentV[i] = srnd.randomFloat();
	}

	for ( i = 0; i < MIXBUFFER_SAMPLES*6; i++ ) {
		origMixBuffer[i] = srnd.randomFloat();
		samples[i] = srnd.randomInt( (1<<16) ) - (1<<15);
	}

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		for ( j = 0; j < MIXBUFFER_SAMPLES*6; j++ ) {
			mixBuffer1[j] = origMixBuffer[j];
		}
		StartRecordTimeLocal( start );
		p_generic-> mixSoundTwoSpeakerMono( mixBuffer1, samples, MIXBUFFER_SAMPLES, lastV, currentV );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic-> mixSoundTwoSpeakerMono()", MIXBUFFER_SAMPLES, bestClocksGeneric );


	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		for ( j = 0; j < MIXBUFFER_SAMPLES*6; j++ ) {
			mixBuffer2[j] = origMixBuffer[j];
		}
		StartRecordTimeLocal( start );
		p_simd-> mixSoundTwoSpeakerMono( mixBuffer2, samples, MIXBUFFER_SAMPLES, lastV, currentV );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < MIXBUFFER_SAMPLES*6; i++ ) {
		if ( CMath::fabs( mixBuffer1[i] - mixBuffer2[i] ) > SOUND_MIX_EPSILON ) {
			break;
		}
	}
	result = ( i >= MIXBUFFER_SAMPLES*6 ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd-> mixSoundTwoSpeakerMono() %s", result ), MIXBUFFER_SAMPLES, bestClocksSIMD, bestClocksGeneric );

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		for ( j = 0; j < MIXBUFFER_SAMPLES*6; j++ ) {
			mixBuffer1[j] = origMixBuffer[j];
		}
		StartRecordTimeLocal( start );
		p_generic-> mixSoundTwoSpeakerStereo( mixBuffer1, samples, MIXBUFFER_SAMPLES, lastV, currentV );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic-> mixSoundTwoSpeakerStereo()", MIXBUFFER_SAMPLES, bestClocksGeneric );


	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		for ( j = 0; j < MIXBUFFER_SAMPLES*6; j++ ) {
			mixBuffer2[j] = origMixBuffer[j];
		}
		StartRecordTimeLocal( start );
		p_simd-> mixSoundTwoSpeakerStereo( mixBuffer2, samples, MIXBUFFER_SAMPLES, lastV, currentV );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < MIXBUFFER_SAMPLES*6; i++ ) {
		if ( CMath::fabs( mixBuffer1[i] - mixBuffer2[i] ) > SOUND_MIX_EPSILON ) {
			break;
		}
	}
	result = ( i >= MIXBUFFER_SAMPLES*6 ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd-> mixSoundTwoSpeakerStereo() %s", result ), MIXBUFFER_SAMPLES, bestClocksSIMD, bestClocksGeneric );


	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		for ( j = 0; j < MIXBUFFER_SAMPLES*6; j++ ) {
			mixBuffer1[j] = origMixBuffer[j];
		}
		StartRecordTimeLocal( start );
		p_generic-> mixSoundSixSpeakerMono( mixBuffer1, samples, MIXBUFFER_SAMPLES, lastV, currentV );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic-> mixSoundSixSpeakerMono()", MIXBUFFER_SAMPLES, bestClocksGeneric );


	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		for ( j = 0; j < MIXBUFFER_SAMPLES*6; j++ ) {
			mixBuffer2[j] = origMixBuffer[j];
		}
		StartRecordTimeLocal( start );
		p_simd-> mixSoundSixSpeakerMono( mixBuffer2, samples, MIXBUFFER_SAMPLES, lastV, currentV );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < MIXBUFFER_SAMPLES*6; i++ ) {
		if ( CMath::fabs( mixBuffer1[i] - mixBuffer2[i] ) > SOUND_MIX_EPSILON ) {
			break;
		}
	}
	result = ( i >= MIXBUFFER_SAMPLES*6 ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd-> mixSoundSixSpeakerMono() %s", result ), MIXBUFFER_SAMPLES, bestClocksSIMD, bestClocksGeneric );

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		for ( j = 0; j < MIXBUFFER_SAMPLES*6; j++ ) {
			mixBuffer1[j] = origMixBuffer[j];
		}
		StartRecordTimeLocal( start );
		p_generic-> mixSoundSixSpeakerStereo( mixBuffer1, samples, MIXBUFFER_SAMPLES, lastV, currentV );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic-> mixSoundSixSpeakerStereo()", MIXBUFFER_SAMPLES, bestClocksGeneric );


	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		for ( j = 0; j < MIXBUFFER_SAMPLES*6; j++ ) {
			mixBuffer2[j] = origMixBuffer[j];
		}
		StartRecordTimeLocal( start );
		p_simd-> mixSoundSixSpeakerStereo( mixBuffer2, samples, MIXBUFFER_SAMPLES, lastV, currentV );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < MIXBUFFER_SAMPLES*6; i++ ) {
		if ( CMath::fabs( mixBuffer1[i] - mixBuffer2[i] ) > SOUND_MIX_EPSILON ) {
			break;
		}
	}
	result = ( i >= MIXBUFFER_SAMPLES*6 ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd-> mixSoundSixSpeakerStereo() %s", result ), MIXBUFFER_SAMPLES, bestClocksSIMD, bestClocksGeneric );


	for ( i = 0; i < MIXBUFFER_SAMPLES*6; i++ ) {
		origMixBuffer[i] = srnd.randomInt( (1<<17) ) - (1<<16);
	}

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		for ( j = 0; j < MIXBUFFER_SAMPLES*6; j++ ) {
			mixBuffer1[j] = origMixBuffer[j];
		}
		StartRecordTimeLocal( start );
		p_generic-> mixedSoundToSamples( outSamples1, mixBuffer1, MIXBUFFER_SAMPLES*6 );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic-> mixedSoundToSamples()", MIXBUFFER_SAMPLES, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		for ( j = 0; j < MIXBUFFER_SAMPLES*6; j++ ) {
			mixBuffer2[j] = origMixBuffer[j];
		}
		StartRecordTimeLocal( start );
		p_simd-> mixedSoundToSamples( outSamples2, mixBuffer2, MIXBUFFER_SAMPLES*6 );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < MIXBUFFER_SAMPLES*6; i++ ) {
		if ( outSamples1[i] != outSamples2[i] ) {
			break;
		}
	}
	result = ( i >= MIXBUFFER_SAMPLES*6 ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd-> mixedSoundToSamples() %s", result ), MIXBUFFER_SAMPLES, bestClocksSIMD, bestClocksGeneric );
}

/*
============
TestMath
============
*/
void TestMath() {
	int i;
	TIME_TYPE start, end, bestClocks;

	Debug::debug(Debug::math,__FUNCTION__) <<"====================================\n" <<endl;

	float tst = -1.0f;
	float tst2 = 1.0f;
	float testvar = 1.0f;
	CRandom rnd;

	bestClocks = 0;
	tst = rnd.randomFloat();
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		tst = fabs( tst );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
		testvar = ( testvar + tst ) * tst;
		tst = rnd.randomFloat();
	}
	CProcClock::printClocks( "            fabs( tst )", 1, bestClocks );

	bestClocks = 0;
	tst = rnd.randomFloat();
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		int tmp = * ( int * ) &tst;
		tmp &= 0x7FFFFFFF;
		tst = * ( float * ) &tmp;
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
		testvar = ( testvar + tst ) * tst;
		tst = rnd.randomFloat();
	}
	CProcClock::printClocks( "    CMath::fabs( tst )", 1, bestClocks );

	bestClocks = 0;
	tst = 10.0f + 100.0f * rnd.randomFloat();
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		tst = sqrt( tst );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
		testvar = ( testvar + tst ) * tst * 0.01f;
		tst = 10.0f + 100.0f * rnd.randomFloat();
	}
	CProcClock::printClocks( "            sqrt( tst )", 1, bestClocks );

	bestClocks = 0;
	tst = rnd.randomFloat();
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		tst = CMath::sqrt( tst );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
		testvar = ( testvar + tst ) * tst;
		tst = rnd.randomFloat();
	}
	CProcClock::printClocks( "    CMath::sqrt( tst )", 1, bestClocks );

	bestClocks = 0;
	tst = rnd.randomFloat();
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		tst = CMath::sqrt16( tst );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
		testvar = ( testvar + tst ) * tst;
		tst = rnd.randomFloat();
	}
	CProcClock::printClocks( "  CMath::sqrt16( tst )", 1, bestClocks );

	bestClocks = 0;
	tst = rnd.randomFloat();
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		tst = CMath::sqrt64( tst );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
		testvar = ( testvar + tst ) * tst;
		tst = rnd.randomFloat();
	}
	CProcClock::printClocks( "  CMath::sqrt64( tst )", 1, bestClocks );

	bestClocks = 0;
	tst = rnd.randomFloat();
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		tst = tst * CMath::rSqrt( tst );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
		testvar = ( testvar + tst ) * tst;
		tst = rnd.randomFloat();
	}
	CProcClock::printClocks( "   CMath::rSqrt( tst )", 1, bestClocks );

	bestClocks = 0;
	tst = rnd.randomFloat();
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		tst = CMath::sin( tst );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
		testvar = ( testvar + tst ) * tst;
		tst = rnd.randomFloat();
	}
	CProcClock::printClocks( "     CMath::sin( tst )", 1, bestClocks );

	bestClocks = 0;
	tst = rnd.randomFloat();
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		tst = CMath::sin16( tst );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
		testvar = ( testvar + tst ) * tst;
		tst = rnd.randomFloat();
	}
	CProcClock::printClocks( "   CMath::sin16( tst )", 1, bestClocks );

	bestClocks = 0;
	tst = rnd.randomFloat();
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		tst = CMath::cos( tst );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
		testvar = ( testvar + tst ) * tst;
		tst = rnd.randomFloat();
	}
	CProcClock::printClocks( "     CMath::cos( tst )", 1, bestClocks );

	bestClocks = 0;
	tst = rnd.randomFloat();
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		tst = CMath::cos16( tst );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
		testvar = ( testvar + tst ) * tst;
		tst = rnd.randomFloat();
	}
	CProcClock::printClocks( "   CMath::cos16( tst )", 1, bestClocks );

	bestClocks = 0;
	tst = rnd.randomFloat();
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		CMath::sincos( tst, tst, tst2 );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
		testvar = ( testvar + tst ) * tst;
		tst = rnd.randomFloat();
	}
	CProcClock::printClocks( "  CMath::sincos( tst )", 1, bestClocks );

	bestClocks = 0;
	tst = rnd.randomFloat();
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		CMath::sincos16( tst, tst, tst2 );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
		testvar = ( testvar + tst ) * tst;
		tst = rnd.randomFloat();
	}
	CProcClock::printClocks( "CMath::sincos16( tst )", 1, bestClocks );

	bestClocks = 0;
	tst = rnd.randomFloat();
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		tst = CMath::tan( tst );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
		testvar = ( testvar + tst ) * tst;
		tst = rnd.randomFloat();
	}
	CProcClock::printClocks( "     CMath::tan( tst )", 1, bestClocks );

	bestClocks = 0;
	tst = rnd.randomFloat();
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		tst = CMath::tan16( tst );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
		testvar = ( testvar + tst ) * tst;
		tst = rnd.randomFloat();
	}
	CProcClock::printClocks( "   CMath::tan16( tst )", 1, bestClocks );

	bestClocks = 0;
	tst = rnd.randomFloat();
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		tst = CMath::asin( tst );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
		testvar = ( testvar + tst ) * tst * ( 1.0f / CMath::PI );
		tst = rnd.randomFloat();
	}
	CProcClock::printClocks( "    CMath::asin( tst )", 1, bestClocks );

	bestClocks = 0;
	tst = rnd.randomFloat();
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		tst = CMath::asin16( tst );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
		testvar = ( testvar + tst ) * tst * ( 1.0f / CMath::PI );
		tst = rnd.randomFloat();
	}
	CProcClock::printClocks( "  CMath::asin16( tst )", 1, bestClocks );

	bestClocks = 0;
	tst = rnd.randomFloat();
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		tst = CMath::acos( tst );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
		testvar = ( testvar + tst ) * tst * ( 1.0f / CMath::PI );
		tst = rnd.randomFloat();
	}
	CProcClock::printClocks( "    CMath::acos( tst )", 1, bestClocks );

	bestClocks = 0;
	tst = rnd.randomFloat();
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		tst = CMath::acos16( tst );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
		testvar = ( testvar + tst ) * tst * ( 1.0f / CMath::PI );
		tst = rnd.randomFloat();
	}
	CProcClock::printClocks( "  CMath::acos16( tst )", 1, bestClocks );

	bestClocks = 0;
	tst = rnd.randomFloat();
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		tst = CMath::atan( tst );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
		testvar = ( testvar + tst ) * tst;
		tst = rnd.randomFloat();
	}
	CProcClock::printClocks( "    CMath::atan( tst )", 1, bestClocks );

	bestClocks = 0;
	tst = rnd.randomFloat();
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		tst = CMath::atan16( tst );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
		testvar = ( testvar + tst ) * tst;
		tst = rnd.randomFloat();
	}
	CProcClock::printClocks( "  CMath::atan16( tst )", 1, bestClocks );

	bestClocks = 0;
	tst = rnd.randomFloat();
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		tst = CMath::pow( 2.7f, tst );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
		testvar = ( testvar + tst ) * tst * 0.1f;
		tst = rnd.randomFloat();
	}
	CProcClock::printClocks( "    CMath::pow( tst )", 1, bestClocks );

	bestClocks = 0;
	tst = rnd.randomFloat();
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		tst = CMath::pow16( 2.7f, tst );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
		testvar = ( testvar + tst ) * tst * 0.1f;
		tst = rnd.randomFloat();
	}
	CProcClock::printClocks( "  CMath::pow16( tst )", 1, bestClocks );

	bestClocks = 0;
	tst = rnd.randomFloat();
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		tst = CMath::exp( tst );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
		testvar = ( testvar + tst ) * tst * 0.1f;
		tst = rnd.randomFloat();
	}
	CProcClock::printClocks( "    CMath::exp( tst )", 1, bestClocks );

	bestClocks = 0;
	tst = rnd.randomFloat();
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		tst = CMath::exp16( tst );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
		testvar = ( testvar + tst ) * tst * 0.1f;
		tst = rnd.randomFloat();
	}
	CProcClock::printClocks( "  CMath::exp16( tst )", 1, bestClocks );

	bestClocks = 0;
	tst = rnd.randomFloat();
	for ( i = 0; i < NUMTESTS; i++ ) {
		tst = fabs( tst ) + 1.0f;
		StartRecordTimeLocal( start );
		tst = CMath::log( tst );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
		testvar = ( testvar + tst ) * tst;
		tst = rnd.randomFloat();
	}
	CProcClock::printClocks( "    CMath::log( tst )", 1, bestClocks );

	bestClocks = 0;
	tst = rnd.randomFloat();
	for ( i = 0; i < NUMTESTS; i++ ) {
		tst = fabs( tst ) + 1.0f;
		StartRecordTimeLocal( start );
		tst = CMath::log16( tst );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
		testvar = ( testvar + tst ) * tst;
		tst = rnd.randomFloat();
	}
	CProcClock::printClocks( "  CMath::log16( tst )", 1, bestClocks );

	Debug::debug(Debug::math,__FUNCTION__) << "testvar = "<< testvar<<endl;

	CMat3D resultMat3;
	CQuaternion fromQuat, toQuat, resultQuat;
	CCompQuaternion cq;
	CEulerAngles ang;

	fromQuat = CEulerAngles( 30, 45, 0 ).toQuat();
	toQuat = CEulerAngles( 45, 0, 0 ).toQuat();
	cq = CEulerAngles( 30, 45, 0 ).toQuat().toQuat();
	ang = CEulerAngles( 30, 40, 50 );

	bestClocks = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		resultMat3 = fromQuat.toMat3();
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
	}
	CProcClock::printClocks( "       CQuaternion::toMat3()", 1, bestClocks );

	bestClocks = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		resultQuat.slerp( fromQuat, toQuat, 0.3f );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
	}
	CProcClock::printClocks( "        CQuaternion::slerp()", 1, bestClocks );

	bestClocks = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		resultQuat = cq.toQuat();
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
	}
	CProcClock::printClocks( "      CCompQuaternion::toQuat()", 1, bestClocks );

	bestClocks = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		resultQuat = ang.toQuat();
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
	}
	CProcClock::printClocks( "     CEulerAngles::toQuat()", 1, bestClocks );

	bestClocks = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {
		StartRecordTimeLocal( start );
		resultMat3 = ang.toMat3();
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocks );
	}
	CProcClock::printClocks( "     CEulerAngles::toMat3()", 1, bestClocks );
}

/*
============
TestNegate
============
*/

// this wasn't previously in the test
void TestNegate() {
	int i;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( float fsrc0[COUNT] );
	ALIGN16( float fsrc1[COUNT] );
	ALIGN16( float fsrc2[COUNT] );

	const char *result;

	CRandom srnd( RANDOM_SEED );

	for ( i = 0; i < COUNT; i++ ) {
		fsrc0[i] = fsrc1[i] = fsrc2[i] = srnd.randomFloat() * 10.0f;
		//fsrc1[i] = srnd.randomFloat() * 10.0f;
	}

	Debug::debug(Debug::math,__FUNCTION__) <<"====================================\n" <<endl;

	bestClocksGeneric = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {

		memcpy( &fsrc1[0], &fsrc0[0], COUNT * sizeof(float) );

		StartRecordTimeLocal( start );
		p_generic->negate16( fsrc1, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->negate16( float[] )", COUNT, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < NUMTESTS; i++ ) {

		memcpy( &fsrc2[0], &fsrc0[0], COUNT * sizeof(float) );

		StartRecordTimeLocal( start );
		p_simd->negate16( fsrc2, COUNT );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < COUNT; i++ ) {
		if ( fsrc1[i] != fsrc2[i] ) {
			break;
		}
	}
	result = ( i >= COUNT ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->negate16( float[] ) %s", result ), COUNT, bestClocksSIMD, bestClocksGeneric );
}


/*
============
CSIMD::Test_f
============
*/
void CSIMD::Test_f( const Util::CCMDLineArgs &args ) {

#ifdef _WIN32
	SetThreadPriority( GetCurrentThread(), THREAD_PRIORITY_TIME_CRITICAL );
#endif /* _WIN32 */

	p_simd = processor;
	p_generic = generic;

	if ( CMyString::length( args.Argv( 1 ) ) != 0 ) {
		cpuid_t cpuid = System::getProcessorId();
		CMyString argString = args.Args();

		argString.replace( " ", "" );

		if ( CMyString::compareInsen( argString, "MMX" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) ) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX\n" <<endl;
				return;
			}
			p_simd = new CSIMD_MMX;
		#if !defined (__GNUC__)
		} else if ( CMyString::compareInsen( argString, "3DNow" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) || !( cpuid & CPUID_3DNOW ) ) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX & 3DNow\n"  <<endl;
				return;
			}
			p_simd = new CSIMD_3DNow;
		#endif
		} else if ( CMyString::compareInsen( argString, "SSE" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) || !( cpuid & CPUID_SSE ) ) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX & SSE\n" <<endl;
				return;
			}
			p_simd = new CSIMD_SSE;
		} else if ( CMyString::compareInsen( argString, "SSE2" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) || !( cpuid & CPUID_SSE ) || !( cpuid & CPUID_SSE2 ) ) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX & SSE & SSE2\n"  <<endl;
				return;
			}
			p_simd = new CSIMD_SSE2;
		} else if ( CMyString::compareInsen( argString, "SSE3" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) || !( cpuid & CPUID_SSE ) || !( cpuid & CPUID_SSE2 ) || !( cpuid & CPUID_SSE3 ) ) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX & SSE & SSE2 & SSE3\n"  <<endl;
				return;
			}
			p_simd = new CSIMD_SSE3();
		}else if ( CMyString::compareInsen( argString, "SSE41" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) || !( cpuid & CPUID_SSE ) || !( cpuid & CPUID_SSE2 ) || !( cpuid & CPUID_SSE3 )|| !( cpuid & CPUID_SSE41 ) ) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX & SSE & SSE2 & SSE3 & SSE41\n"  <<endl;
				return;
			}
			p_simd = new CSIMD_SSE41;//CCPUID_SSE41();
		}else if ( CMyString::compareInsen( argString, "SSE42" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) || !( cpuid & CPUID_SSE ) || !( cpuid & CPUID_SSE2 ) || !( cpuid & CPUID_SSE3 )|| !( cpuid & CPUID_SSE41 ) || !( cpuid & CPUID_SSE42 )) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX & SSE & SSE2 & SSE3 & SSE41 & SSE41\n"  <<endl;
				return;
			}
			p_simd = new CSIMD_SSE42;//CSIMD_SSE42();
		}else if ( CMyString::compareInsen( argString, "AVX" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) || !( cpuid & CPUID_SSE ) || !( cpuid & CPUID_SSE2 ) || !( cpuid & CPUID_SSE3 )|| !( cpuid & CPUID_SSE41 ) || !( cpuid & CPUID_SSE42 ) || !( cpuid & CPUID_AVX )) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX & SSE & SSE2 & SSE3 & SSE41 & SSE41 & AVX\n"  <<endl;
				return;
			}
			p_simd = new CSIMD_AVX;
		} /*else if ( CMyString::compareInsen( argString, "AltiVec" ) == 0 ) {
			if ( !( cpuid & CPUID_ALTIVEC ) ) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support AltiVec\n"  <<endl;
				return;
			}
			p_simd = new CSIMD_AltiVec();
		}*/ else {
			Debug::debug(Debug::math,__FUNCTION__) << "invalid argument, use: MMX, 3DNow, SSE, SSE2, SSE3, AltiVec\n"  <<endl;
			return;
		}
	}

	//idLib::common->SetRefreshOnPrint( true );

	Debug::debug(Debug::math,__FUNCTION__) << "using %s for SIMD processing: "<< p_simd->getName() <<endl;

	CProcClock::getBaseClocks();

	TestMath();
	TestAdd();
	TestSub();
	TestMul();
	TestDiv();
	TestMulAdd();
	TestMulSub();
	TestDot();
	TestCompare();
	TestMinMax();
	TestClamp();
	TestMemcpy();
	TestMemset();
	TestNegate();

	TestMatXMultiplyVecX();
	TestMatXMultiplyAddVecX();
	TestMatXTransposeMultiplyVecX();
	TestMatXTransposeMultiplyAddVecX();
	TestMatXMultiplyMatX();
	TestMatXTransposeMultiplyMatX();
	TestMatXLowerTriangularsolve();
	TestMatXLowerTriangularsolveTranspose();
	TestMatXLDLTFactor();

	Debug::debug(Debug::math,__FUNCTION__) <<"====================================\n" <<endl;

	TestblendJoints();
	TestconvertJointQuatsToJointMats();
	TestconvertJointMatsToJointQuats();
	TesttransformJoints();
	TestuntransformJoints();
	TesttransformVerts();
	TesttracePointCull();
	TestdecalPointCull();
	TestoverlayPointCull();
	TestderiveTriPlanes();
	Test_deriveTangents();
	//s TestDeriveUnsmoothedTangents();
	Test_normalizeTangents();
	TestGetTextureSpaceLightVectors();
	TestGetSpecularTextureCoords();
	Test_createShadowCache();

	Debug::debug(Debug::math,__FUNCTION__) <<"====================================\n" <<endl;

	TestSoundUpSampling();
	TestSoundMixing();

	//idLib::common->SetRefreshOnPrint( false );

	if ( p_simd != processor ) {
		delete p_simd;
	}
	p_simd = NULL;
	p_generic = NULL;

#ifdef _WIN32
	SetThreadPriority( GetCurrentThread(), THREAD_PRIORITY_NORMAL );
#endif /* _WIN32 */
}
//====================TEST2===========================================

CVec3D out, in1, in2;
ALIGNTO16 CVec3D outAligned;
ALIGNTO16 CVec3D in1Aligned;
ALIGNTO16 CVec3D in2Aligned;

static void reset3D()
{
	out.toZero();
	outAligned.toZero();

	in1.x = 10.0f;
	in1.y = 20.0f;
	in1.z = 30.0f;

	in2.x = 1.0f;
	in2.y = 2.0f;
	in2.z = 3.0f;

	in1Aligned.x = 10.0f;
	in1Aligned.y = 20.0f;
	in1Aligned.z = 30.0f;

	in2Aligned.x = 1.0f;
	in2Aligned.y = 2.0f;
	in2Aligned.z = 3.0f;

}



static void TestVector3DSum()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec3D tst,dst;
	const char *result;
	in1.toZero();
	in2.toZero();
	reset3D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector3D Sum()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector3D_Sum(&in1, &in2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = in1;

		CProcClock::printClocks( va( "generic->vector3D_Sum "), 1, bestClocksGeneric );
		reset3D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector3D_Sum(&in1, &in2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = in1;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector3D_Sum  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestVector3DSumOf()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec3D tst,dst;
	const char *result;
	in1.toZero();
	in2.toZero();
	reset3D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector3D SumOf()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector3D_SumOf(&out, &in1, &in2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = out;

		CProcClock::printClocks( va( "generic->vector3D_SumOf "), 1, bestClocksGeneric );
		reset3D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector3D_SumOf(&out, &in1, &in2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = out;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector3D_SumOf  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestVector3DDiff()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec3D tst,dst;
	const char *result;
	in1.toZero();
	in2.toZero();
	reset3D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector3D Diff()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector3D_Diff(&in1, &in2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = in1;

		CProcClock::printClocks( va( "generic->vector3D_Diff "), 1, bestClocksGeneric );
		reset3D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector3D_Diff(&in1, &in2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = in1;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector3D_Diff  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestVector3DDiffOf()
{

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec3D tst,dst;
	const char *result;
	in1.toZero();
	in2.toZero();
	reset3D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector3D DiffOf()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector3D_DiffOf(&out, &in1, &in2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = out;

		CProcClock::printClocks( va( "generic->vector3D_DiffOf "), 1, bestClocksGeneric );
		reset3D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector3D_DiffOf(&out, &in1, &in2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = out;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->Vector3DD_DiffOf  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}
static void TestVector3DScale()
{

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec3D tst,dst;
	const char *result;
	in1.toZero();
	in2.toZero();
	reset3D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector3D scale()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector3D_Scale(&in1, 100.0f);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = in1;

		CProcClock::printClocks( va( "generic->vector3D_Scale "), 1, bestClocksGeneric );
		reset3D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector3D_Scale(&in1, 100.0f);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = in1;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector3D_Scale  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestVector3DScaleOf()
{

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec3D tst,dst;
	const char *result;
	in1.toZero();
	in2.toZero();
	reset3D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector3D ScaleOf()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector3D_ScaleOf(&out, &in1, 100.0f);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = out;

		CProcClock::printClocks( va( "generic->vector3D_ScaleOf "), 1, bestClocksGeneric );
		reset3D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector3D_ScaleOf(&out, &in1, 100.0f);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = out;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector3D_ScaleOf  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}


static void TestVector3DDot()
{

	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	CVec3D tst,dst;
	const char *result;
	in1.toZero();
	in2.toZero();
	reset3D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector3D dot()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector3D_Dot(&in1, &in2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = in1;

		CProcClock::printClocks( va( "generic->vector3D_Dot "), 1, bestClocksGeneric );
		reset3D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector3D_Dot(&in1, &in2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = in1;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector3D_Dot  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );




}


static void TestCVec3DLengthSq()
{

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	const char *result;
	float resultado1,resultado2;
	in1.toZero();
	in2.toZero();
	reset3D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector3D getLengthSqr()..."<<endl;


			StartRecordTimeLocal( start );
			resultado1 = p_generic->vector3D_LengthSq(&in2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		CProcClock::printClocks( va( "generic->vector3D_LengthSq "), 1, bestClocksGeneric );
		reset3D();
		bestClocksSIMD = 0;


			StartRecordTimeLocal( start );
			printf("Before Teste\n");

			resultado2 = p_simd->vector3D_LengthSq(&in2);
			printf("After Teste\n");
            StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );


		result = resultado1 == resultado2 ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector3D_LengthSq  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );

}

static void TestCVec3DLength()
{
	reset3D();

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	const char *result;
	float resultado1,resultado2;
	in1.toZero();
	in2.toZero();
	reset3D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector3D getLenght()..."<<endl;


			StartRecordTimeLocal( start );
			resultado1 = p_generic->vector3D_Length(&in2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );


		CProcClock::printClocks( va( "generic->vector3D_Length "), 1, bestClocksGeneric );
		reset3D();
		bestClocksSIMD = 0;


			StartRecordTimeLocal( start );
			resultado2=p_simd->vector3D_Length(&in2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );


		result = resultado1 == resultado2 ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector3D_Length  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );

}

static void TestVector3DCross()
{
	in1.toZero();
	in2.toZero();

	reset3D();
	in1.x = 1.0f; in1.y = 0.0f; in1.z = 0.0f;
	in2.x = 0.0f; in2.y = 0.0f; in2.z = 1.0f;


	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec3D tst,dst;
	const char *result;

	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector3D cross()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector3D_Cross(&in1, &in2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = in1;

		CProcClock::printClocks( va( "generic->vector3D_Cross "), 1, bestClocksGeneric );
		reset3D();
		in1.x = 1.0f; in1.y = 0.0f; in1.z = 0.0f;
		in2.x = 0.0f; in2.y = 0.0f; in2.z = 1.0f;
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector3D_Cross(&in1, &in2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = in1;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector3D_Cross  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );





}

static void TestVector3DCrossOf()
{
	reset3D();
	in1.toZero();
	in2.toZero();
	in1.x = 1.0f; in1.y = 0.0f; in1.z = 0.0f;
	in2.x = 0.0f; in2.y = 0.0f; in2.z = 1.0f;

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec3D tst,dst;
	const char *result;


	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector3D CrossOf()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector3D_CrossOf(&out, &in1, &in2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = out;

		CProcClock::printClocks( va( "generic->vector3D_CrossOf "), 1, bestClocksGeneric );
		reset3D();

		in1.x = 1.0f; in1.y = 0.0f; in1.z = 0.0f;
		in2.x = 0.0f; in2.y = 0.0f; in2.z = 1.0f;

		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector3D_CrossOf(&out, &in1, &in2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = out;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector3D_CrossOf  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestCVec3DNormalize()
{

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec3D tst,dst;
	const char *result;
	in1.toZero();
	in2.toZero();
	reset3D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector3D  toNormal()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector3D_Normalize(&in1);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = in1;

		CProcClock::printClocks( va( "generic->vector3D_Normalize "), 1, bestClocksGeneric );
		reset3D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector3D_Normalize(&in1);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = in1;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector3D_Normalize  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



	;
}

static void TestCVec3DNormalizeOf()
{

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec3D tst,dst;
	const char *result;
	in1.toZero();
	in2.toZero();
	reset3D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector3D NormalizeOf()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector3D_NormalizeOf(&out, &in1);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = out;

		CProcClock::printClocks( va( "generic->vector3D_NormalizeOf "), 1, bestClocksGeneric );
		reset3D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector3D_NormalizeOf(&out, &in1);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = out;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector3D_NormalizeOf  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestCVec3DDistance()
{
	reset3D();

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	const char *result;
	float resultado1,resultado2;
	in1.toZero();
	in2.toZero();
	reset3D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector3D distance()..."<<endl;


			StartRecordTimeLocal( start );
			resultado1 = p_generic->vector3D_Distance(&in1,&in2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );


		CProcClock::printClocks( va( "generic->vector3D_Distance "), 1, bestClocksGeneric );
		reset3D();
		bestClocksSIMD = 0;


			StartRecordTimeLocal( start );
			resultado2=p_simd->vector3D_Distance(&in1,&in2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );


		result = resultado1 == resultado2 ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector3D_Distance  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );

}



static void TestVector3DAlignedSum()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec3D tst,dst;
	const char *result;
	in1Aligned.toZero();
	in2Aligned.toZero();
	reset3D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector3D Aligned Sum()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector3D_AlignedSum(&in1Aligned, &in2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = in1Aligned;

		CProcClock::printClocks( va( "generic->vector3D_AlignedSum "), 1, bestClocksGeneric );
		reset3D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector3D_AlignedSum(&in1Aligned, &in2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = in1Aligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector3D_AlignedSum  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestVector3DAlignedSumOf()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec3D tst,dst;
	const char *result;
	in1Aligned.toZero();
	in2Aligned.toZero();
	reset3D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector3D Aligned SumOf()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector3D_AlignedSumOf(&outAligned, &in1Aligned, &in2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = outAligned;

		CProcClock::printClocks( va( "generic->vector3D_AlignedSumOf "), 1, bestClocksGeneric );
		reset3D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector3D_AlignedSumOf(&outAligned, &in1Aligned, &in2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = outAligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector3D_AlignedSumOf  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestVector3DAlignedDiff()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec3D tst,dst;
	const char *result;
	in1Aligned.toZero();
	in2Aligned.toZero();
	reset3D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector3D Aligned Diff()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector3D_AlignedDiff(&in1Aligned, &in2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = in1Aligned;

		CProcClock::printClocks( va( "generic->vector3D_AlignedDiff "), 1, bestClocksGeneric );
		reset3D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector3D_AlignedDiff(&in1Aligned, &in2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = in1Aligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector3D_AlignedDiff  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestVector3DAlignedDiffOf()
{

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec3D tst,dst;
	const char *result;
	in1Aligned.toZero();
	in2Aligned.toZero();
	reset3D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector3D Aligned DiffOf()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector3D_AlignedDiffOf(&outAligned, &in1Aligned, &in2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = outAligned;

		CProcClock::printClocks( va( "generic->vector3D_AlignedDiffOf "), 1, bestClocksGeneric );
		reset3D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector3D_AlignedDiffOf(&outAligned, &in1Aligned, &in2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = outAligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector3D_AlignedDiffOf  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}
static void TestVector3DAlignedScale()
{

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec3D tst,dst;
	const char *result;
	in1Aligned.toZero();
	in2Aligned.toZero();
	reset3D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector3D Aligned scale()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector3D_AlignedScale(&in1Aligned, 100.0f);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = in1Aligned;

		CProcClock::printClocks( va( "generic->vector3D_AlignedScale "), 1, bestClocksGeneric );
		reset3D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector3D_AlignedScale(&in1Aligned, 100.0f);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = in1Aligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector3D_AlignedScale  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestVector3DAlignedScaleOf()
{

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec3D tst,dst;
	const char *result;
	in1Aligned.toZero();
	in2Aligned.toZero();
	reset3D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector3D Aligned ScaleOf()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector3D_AlignedScaleOf(&outAligned, &in1Aligned, 100.0f);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = outAligned;

		CProcClock::printClocks( va( "generic->vector3D_AlignedScaleOf "), 1, bestClocksGeneric );
		reset3D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector3D_AlignedScaleOf(&outAligned, &in1Aligned, 100.0f);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = outAligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector3D_AlignedScaleOf  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}


static void TestVector3DAlignedDot()
{

	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	CVec3D tst,dst;
	const char *result;
	in1Aligned.toZero();
	in2Aligned.toZero();
	reset3D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector3D Aligned dot()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector3D_AlignedDot(&in1Aligned, &in2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = in1Aligned;

		CProcClock::printClocks( va( "generic->vector3D_AlignedDot "), 1, bestClocksGeneric );
		reset3D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector3D_AlignedDot(&in1Aligned, &in2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = in1Aligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector3D_AlignedDot  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );




}


static void TestVector3DAlignedLengthSq()
{

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	const char *result;
	float resultado1,resultado2;
	in1Aligned.toZero();
	in2Aligned.toZero();
	reset3D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector3D Aligned getLengthSqr()..."<<endl;


			StartRecordTimeLocal( start );
			resultado1 = p_generic->vector3D_AlignedLengthSq(&in2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		CProcClock::printClocks( va( "generic->vector3D_AlignedLengthSq "), 1, bestClocksGeneric );
		reset3D();
		bestClocksSIMD = 0;


			StartRecordTimeLocal( start );
			printf("Before Teste\n");

			resultado2 = p_simd->vector3D_AlignedLengthSq(&in2Aligned);
			printf("After Teste\n");
            StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );


		result = resultado1 == resultado2 ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector3D_AlignedLengthSq  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );

}

static void TestVector3DAlignedLength()
{
	reset3D();

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	const char *result;
	float resultado1,resultado2;
	in1Aligned.toZero();
	in2Aligned.toZero();
	reset3D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector3D Aligned getLenght()..."<<endl;


			StartRecordTimeLocal( start );
			resultado1 = p_generic->vector3D_AlignedLength(&in2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );


		CProcClock::printClocks( va( "generic->vector3D_AlignedLength "), 1, bestClocksGeneric );
		reset3D();
		bestClocksSIMD = 0;


			StartRecordTimeLocal( start );
			resultado2=p_simd->vector3D_AlignedLength(&in2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );


		result = resultado1 == resultado2 ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector3D_AlignedLength  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );

}

static void TestVector3DAlignedCross()
{
	in1Aligned.toZero();
	in2Aligned.toZero();

	reset3D();
	in1Aligned.x = 1.0f; in1Aligned.y = 0.0f; in1Aligned.z = 0.0f;
	in2Aligned.x = 0.0f; in2Aligned.y = 0.0f; in2Aligned.z = 1.0f;


	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec3D tst,dst;
	const char *result;

	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector3D Aligned cross()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector3D_AlignedCross(&in1Aligned, &in2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = in1Aligned;

		CProcClock::printClocks( va( "generic->vector3D_AlignedCross "), 1, bestClocksGeneric );
		reset3D();
		in1Aligned.x = 1.0f; in1Aligned.y = 0.0f; in1Aligned.z = 0.0f;
		in2Aligned.x = 0.0f; in2Aligned.y = 0.0f; in2Aligned.z = 1.0f;
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector3D_AlignedCross(&in1Aligned, &in2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = in1Aligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector3D_AlignedCross  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );





}

static void TestVector3DAlignedCrossOf()
{
	reset3D();
	in1Aligned.toZero();
	in2Aligned.toZero();
	in1Aligned.x = 1.0f; in1Aligned.y = 0.0f; in1Aligned.z = 0.0f;
	in2Aligned.x = 0.0f; in2Aligned.y = 0.0f; in2Aligned.z = 1.0f;

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec3D tst,dst;
	const char *result;


	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector3D Aligned CrossOf()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector3D_AlignedCrossOf(&outAligned, &in1Aligned, &in2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = outAligned;

		CProcClock::printClocks( va( "generic->vector3D_AlignedCrossOf "), 1, bestClocksGeneric );
		reset3D();

		in1Aligned.x = 1.0f; in1Aligned.y = 0.0f; in1Aligned.z = 0.0f;
		in2Aligned.x = 0.0f; in2Aligned.y = 0.0f; in2Aligned.z = 1.0f;

		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector3D_AlignedCrossOf(&outAligned, &in1Aligned, &in2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = outAligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector3D_AlignedCrossOf  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestCVec3DAlignedNormalize()
{

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec3D tst,dst;
	const char *result;
	in1Aligned.toZero();
	in2Aligned.toZero();
	reset3D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector3D Aligned toNormal()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector3D_AlignedNormalize(&in1Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = in1Aligned;

		CProcClock::printClocks( va( "generic->vector3D_AlignedNormalize "), 1, bestClocksGeneric );
		reset3D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector3D_AlignedNormalize(&in1Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = in1Aligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector3D_AlignedNormalize  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



	;
}

static void TestCVec3DAlignedNormalizeOf()
{

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec3D tst,dst;
	const char *result;
	in1Aligned.toZero();
	in2Aligned.toZero();
	reset3D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector3D Aligned NormalizeOf()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector3D_AlignedNormalizeOf(&outAligned, &in1Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = outAligned;

		CProcClock::printClocks( va( "generic->vector3D_AlignedNormalizeOf "), 1, bestClocksGeneric );
		reset3D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector3D_AlignedNormalizeOf(&outAligned, &in1Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = outAligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector3D_AlignedNormalizeOf  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestCVec3DAlignedDistance()
{
	reset3D();

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	const char *result;
	float resultado1,resultado2;
	in1Aligned.toZero();
	in2Aligned.toZero();
	reset3D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector3D Aligned distance()..."<<endl;


			StartRecordTimeLocal( start );
			resultado1 = p_generic->vector3D_AlignedDistance(&in1Aligned,&in2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );


		CProcClock::printClocks( va( "generic->vector3D_AlignedDistance "), 1, bestClocksGeneric );
		reset3D();
		bestClocksSIMD = 0;


			StartRecordTimeLocal( start );
			resultado2=p_simd->vector3D_AlignedDistance(&in1Aligned,&in2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );


		result = resultado1 == resultado2 ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector3D_AlignedDistance  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );

}
















//=======================================================================================
/*
============
CSIMD::Test3D_f
============
*/

void CSIMD::TesteLenght3D(const Util::CCMDLineArgs &args){
	p_simd = processor;
	p_generic = generic;
	if ( CMyString::length( args.Argv( 1 ) ) != 0 ) {
		cpuid_t cpuid = System::getProcessorId();
		CMyString argString = args.Args();

		argString.replace( " ", "" );

		if ( CMyString::compareInsen( argString, "MMX" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) ) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX\n" <<endl;
				return;
			}
			p_simd = new CSIMD_MMX;
			#if !defined (__GNUC__)
		} else if ( CMyString::compareInsen( argString, "3DNow" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) || !( cpuid & CPUID_3DNOW ) ) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX & 3DNow\n"  <<endl;
				return;
			}
			p_simd = new CSIMD_3DNow;
		#endif
		} else if ( CMyString::compareInsen( argString, "SSE" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) || !( cpuid & CPUID_SSE ) ) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX & SSE\n" <<endl;
				return;
			}
			p_simd = new CSIMD_SSE;
		} else if ( CMyString::compareInsen( argString, "SSE2" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) || !( cpuid & CPUID_SSE ) || !( cpuid & CPUID_SSE2 ) ) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX & SSE & SSE2\n"  <<endl;
				return;
			}
			p_simd = new CSIMD_SSE2;
		} else if ( CMyString::compareInsen( argString, "SSE3" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) || !( cpuid & CPUID_SSE ) || !( cpuid & CPUID_SSE2 ) || !( cpuid & CPUID_SSE3 ) ) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX & SSE & SSE2 & SSE3\n"  <<endl;
				return;
			}
			p_simd = new CSIMD_SSE3();
		}else if ( CMyString::compareInsen( argString, "SSE41" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) || !( cpuid & CPUID_SSE ) || !( cpuid & CPUID_SSE2 ) || !( cpuid & CPUID_SSE3 )|| !( cpuid & CPUID_SSE41 ) ) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX & SSE & SSE2 & SSE3 & SSE41\n"  <<endl;
				return;
			}
			p_simd = new CSIMD_SSE41;//CCPUID_SSE41();
		}else if ( CMyString::compareInsen( argString, "SSE42" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) || !( cpuid & CPUID_SSE ) || !( cpuid & CPUID_SSE2 ) || !( cpuid & CPUID_SSE3 )|| !( cpuid & CPUID_SSE41 ) || !( cpuid & CPUID_SSE42 )) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX & SSE & SSE2 & SSE3 & SSE41 & SSE41\n"  <<endl;
				return;
			}
			p_simd = new CSIMD_SSE42;//CSIMD_SSE42();
		}else if ( CMyString::compareInsen( argString, "AVX" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) || !( cpuid & CPUID_SSE ) || !( cpuid & CPUID_SSE2 ) || !( cpuid & CPUID_SSE3 )|| !( cpuid & CPUID_SSE41 ) || !( cpuid & CPUID_SSE42 ) || !( cpuid & CPUID_AVX )) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX & SSE & SSE2 & SSE3 & SSE41 & SSE41 & AVX\n"  <<endl;
				return;
			}
			p_simd = new CSIMD_AVX;
		} /*else if ( CMyString::compareInsen( argString, "AltiVec" ) == 0 ) {
			if ( !( cpuid & CPUID_ALTIVEC ) ) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support AltiVec\n"  <<endl;
				return;
			}
			p_simd = new CSIMD_AltiVec();
		} */else {
			Debug::debug(Debug::math,__FUNCTION__) << "invalid argument, use: MMX, 3DNow, SSE, SSE2, SSE3, AltiVec\n"  <<endl;
			return;
		}
	}
	Debug::debug(Debug::math,__FUNCTION__) << "using %s for SIMD processing: "<< p_simd->getName() <<endl;

	CProcClock::getBaseClocks();
	//testes=============================================

	TestCVec3DLength();

}


void CSIMD::Test3D_f( const Util::CCMDLineArgs &args ){

	#ifdef _WIN32
	SetThreadPriority( GetCurrentThread(), THREAD_PRIORITY_TIME_CRITICAL );
#endif /* _WIN32 */

	p_simd = processor;
	p_generic = generic;

	if ( CMyString::length( args.Argv( 1 ) ) != 0 ) {
		cpuid_t cpuid = System::getProcessorId();
		CMyString argString = args.Args();

		argString.replace( " ", "" );

		if ( CMyString::compareInsen( argString, "MMX" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) ) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX\n" <<endl;
				return;
			}
			p_simd = new CSIMD_MMX;
			#if !defined (__GNUC__)
		} else if ( CMyString::compareInsen( argString, "3DNow" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) || !( cpuid & CPUID_3DNOW ) ) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX & 3DNow\n"  <<endl;
				return;
			}
			p_simd = new CSIMD_3DNow;
		#endif
		} else if ( CMyString::compareInsen( argString, "SSE" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) || !( cpuid & CPUID_SSE ) ) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX & SSE\n" <<endl;
				return;
			}
			p_simd = new CSIMD_SSE;
		} else if ( CMyString::compareInsen( argString, "SSE2" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) || !( cpuid & CPUID_SSE ) || !( cpuid & CPUID_SSE2 ) ) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX & SSE & SSE2\n"  <<endl;
				return;
			}
			p_simd = new CSIMD_SSE2;
		} else if ( CMyString::compareInsen( argString, "SSE3" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) || !( cpuid & CPUID_SSE ) || !( cpuid & CPUID_SSE2 ) || !( cpuid & CPUID_SSE3 ) ) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX & SSE & SSE2 & SSE3\n"  <<endl;
				return;
			}
			p_simd = new CSIMD_SSE3();
		}else if ( CMyString::compareInsen( argString, "SSE41" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) || !( cpuid & CPUID_SSE ) || !( cpuid & CPUID_SSE2 ) || !( cpuid & CPUID_SSE3 )|| !( cpuid & CPUID_SSE41 ) ) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX & SSE & SSE2 & SSE3 & SSE41\n"  <<endl;
				return;
			}
			p_simd = new CSIMD_SSE41;
		}else if ( CMyString::compareInsen( argString, "SSE42" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) || !( cpuid & CPUID_SSE ) || !( cpuid & CPUID_SSE2 ) || !( cpuid & CPUID_SSE3 )|| !( cpuid & CPUID_SSE41 ) || !( cpuid & CPUID_SSE42 )) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX & SSE & SSE2 & SSE3 & SSE41 & SSE41\n"  <<endl;
				return;
			}
			p_simd = new CSIMD_SSE42;
		}else if ( CMyString::compareInsen( argString, "AVX" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) || !( cpuid & CPUID_SSE ) || !( cpuid & CPUID_SSE2 ) || !( cpuid & CPUID_SSE3 )|| !( cpuid & CPUID_SSE41 ) || !( cpuid & CPUID_SSE42 ) || !( cpuid & CPUID_AVX )) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX & SSE & SSE2 & SSE3 & SSE41 & SSE41 & AVX\n"  <<endl;
				return;
			}
			p_simd = new CSIMD_AVX;
		} /*else if ( CMyString::compareInsen( argString, "AltiVec" ) == 0 ) {
			if ( !( cpuid & CPUID_ALTIVEC ) ) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support AltiVec\n"  <<endl;
				return;
			}
			p_simd = new CSIMD_AltiVec();
		} */else {
			Debug::debug(Debug::math,__FUNCTION__) << "invalid argument, use: MMX, 3DNow, SSE, SSE2, SSE3, AltiVec\n"  <<endl;
			return;
		}
	}

	//idLib::common->SetRefreshOnPrint( true );

	Debug::debug(Debug::math,__FUNCTION__) << "using %s for SIMD processing: "<< p_simd->getName() <<endl;

	CProcClock::getBaseClocks();
	//testes=============================================

	TestVector3DSum();
	TestVector3DSumOf();
	TestVector3DDiff();
	TestVector3DDiffOf();
	TestVector3DScale();
	TestVector3DScaleOf();
	TestVector3DDot();
	TestCVec3DLengthSq();
	TestCVec3DLength();
	TestVector3DCross();
	TestVector3DCrossOf();
	TestCVec3DNormalize();
	TestCVec3DNormalizeOf();
	TestCVec3DDistance();

	TestVector3DAlignedSum();
	TestVector3DAlignedSumOf();
	TestVector3DAlignedDiff();
	TestVector3DAlignedDiffOf();
	TestVector3DAlignedScale();
	TestVector3DAlignedScaleOf();
	TestVector3DAlignedDot();
	TestVector3DAlignedLengthSq();
	TestVector3DAlignedLength();
	TestVector3DAlignedCross();
	TestVector3DAlignedCrossOf();
	TestCVec3DAlignedNormalize();
	TestCVec3DAlignedNormalizeOf();
	TestCVec3DAlignedDistance();
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

//============ TEST CVec4D=================

CVec4D out4d, in4d1, in4d2;
ALIGNTO16  CVec4D out4dAligned;
ALIGNTO16  CVec4D in4d1Aligned;
ALIGNTO16  CVec4D in4d2Aligned;
static void reset4D()
{
	out4d.x = out4d.y = out4d.z= out4d.w = 0.0f;

	in4d1.x = 10.0f;
	in4d1.y = 20.0f;
	in4d1.z = 30.0f;
	in4d1.w = 40.0f;

	in4d2.x = 1.0f;
	in4d2.y = 2.0f;
	in4d2.z = 3.0f;
	in4d2.w = 4.0f;

	in4d1Aligned.x = 10.0f;
	in4d1Aligned.y = 20.0f;
	in4d1Aligned.z = 30.0f;
	in4d1Aligned.w = 40.0f;

	in4d2Aligned.x = 1.0f;
	in4d2Aligned.y = 2.0f;
	in4d2Aligned.z = 3.0f;
	in4d2Aligned.w = 4.0f;

}



static void TestVector4DSum()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec4D tst,dst;
	const char *result;
	in4d1.toZero();
	in4d2.toZero();
	reset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector4D Sum()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector4D_Sum(&in4d1, &in4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = in4d1;

		CProcClock::printClocks( va( "generic->vector4D_Sum "), 1, bestClocksGeneric );
		reset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector4D_Sum(&in4d1, &in4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = in4d1;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector4D_Sum  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}



static void TestVector4DSumOf()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec4D tst,dst;
	const char *result;
	in4d1.toZero();
	in4d2.toZero();
	reset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector4D SumOf()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector4D_SumOf(&out4d, &in4d1, &in4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = out4d;

		CProcClock::printClocks( va( "generic->vector4D_SumOf "), 1, bestClocksGeneric );
		reset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector4D_SumOf(&out4d, &in4d1, &in4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = out4d;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector4D_SumOf  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestVector4DDiff()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec4D tst,dst;
	const char *result;
	in4d1.toZero();
	in4d2.toZero();
	reset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector4D Diff()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector4D_Diff(&in4d1, &in4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = in4d1;

		CProcClock::printClocks( va( "generic->vector4D_Diff "), 1, bestClocksGeneric );
		reset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector4D_Diff(&in4d1, &in4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = in4d1;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector4D_Diff  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestVector4DDiffOf()
{

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec4D tst,dst;
	const char *result;
	in4d1.toZero();
	in4d2.toZero();
	reset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector4D DiffOf()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector4D_DiffOf(&out4d, &in4d1, &in4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = out4d;

		CProcClock::printClocks( va( "generic->vector4D_DiffOf "), 1, bestClocksGeneric );
		reset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector4D_DiffOf(&out4d, &in4d1, &in4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = out4d;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector4D_DiffOf  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}
static void TestVector4DScale()
{

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec4D tst,dst;
	const char *result;
	in4d1.toZero();
	in4d2.toZero();
	reset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector4D scale()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector4D_Scale(&in4d1, 100.0f);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = in4d1;

		CProcClock::printClocks( va( "generic->vector4D_Scale "), 1, bestClocksGeneric );
		reset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector4D_Scale(&in4d1, 100.0f);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = in4d1;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector4D_Scale  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestVector4DScaleOf()
{

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec4D tst,dst;
	const char *result;
	in4d1.toZero();
	in4d2.toZero();
	reset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector4D ScaleOf()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector4D_ScaleOf(&out4d, &in4d1, 100.0f);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = out4d;

		CProcClock::printClocks( va( "generic->vector4D_ScaleOf "), 1, bestClocksGeneric );
		reset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector4D_ScaleOf(&out4d, &in4d1, 100.0f);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = out4d;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector4D_ScaleOf  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}



static void TestVector4DDot()
{
	CVec4D in4d1( 1.0f, 2.0f, 3.0f, 4.0f );
	CVec4D in4d2( 5.0f, 6.0f, 7.0f, 8.0f );

	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	CVec4D tst,dst;
	const char *result;
	in4d1.toZero();
	in4d2.toZero();
	reset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector4D Sum()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic-> vector4D_Dot(&in4d1, &in4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = in4d1;

		CProcClock::printClocks( va( "generic->vector4D_Dot "), 1, bestClocksGeneric );
		reset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector4D_Dot(&in4d1, &in4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = in4d1;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector4D_Dot  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );


	/* Should be 5 + 12 + 21 + 32 = 70 */


}


static void TestCVec4DLengthSq()
{

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	const char *result;
	float resultado1,resultado2;
	in4d1.toZero();
	in4d2.toZero();
	reset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector4D getLengthSqr()..."<<endl;


			StartRecordTimeLocal( start );
			resultado1 = p_generic->vector4D_LengthSq(&in4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		CProcClock::printClocks( va( "generic->vector4D_LengthSq "), 1, bestClocksGeneric );
		reset4D();
		bestClocksSIMD = 0;


			StartRecordTimeLocal( start );
			resultado2 = p_simd->vector4D_LengthSq(&in4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );


		result = resultado1 == resultado2 ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector4D_LengthSq  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );

}

static void TestCVec4DLength()
{
	reset4D();

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	const char *result;
	float resultado1,resultado2;
	in4d1.toZero();
	in4d2.toZero();
	reset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector4D getLenght()..."<<endl;


			StartRecordTimeLocal( start );
			resultado1 = p_generic->vector4D_Length(&in4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );


		CProcClock::printClocks( va( "generic->vector4D_Length "), 1, bestClocksGeneric );
		reset4D();
		bestClocksSIMD = 0;


			StartRecordTimeLocal( start );
			resultado2=p_simd->vector4D_Length(&in4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );


		result = resultado1 == resultado2 ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector4D_Length  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );

}


static void TestCVec4DNormalize()
{

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec4D tst,dst;
	const char *result;
	in4d1.toZero();
	in4d2.toZero();
	reset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector4D  toNormal()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector4D_Normalize(&in4d1);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = in4d1;

		CProcClock::printClocks( va( "generic->vector4D_Normalize "), 1, bestClocksGeneric );
		reset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector4D_Normalize(&in4d1);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = in4d1;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector4D_Normalize  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



	;
}

static void TestCVec4DNormalizeOf()
{

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CVec4D tst,dst;
	const char *result;
	in4d1.toZero();
	in4d2.toZero();
	reset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector4D NormalizeOf()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector4D_NormalizeOf(&out4d, &in4d1);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = out4d;

		CProcClock::printClocks( va( "generic->vector4D_NormalizeOf "), 1, bestClocksGeneric );
		reset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector4D_NormalizeOf(&out4d, &in4d1);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = out4d;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector4D_NormalizeOf  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestCVec4DDistance()
{
	reset4D();

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	const char *result;
	float resultado1,resultado2;
	in4d1.toZero();
	in4d2.toZero();
	reset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector4D distance()..."<<endl;


			StartRecordTimeLocal( start );
			resultado1 = p_generic->vector4D_Distance(&in4d1,&in4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );


		CProcClock::printClocks( va( "generic->vector4D_Distance "), 1, bestClocksGeneric );
		reset4D();
		bestClocksSIMD = 0;


			StartRecordTimeLocal( start );
			resultado2=p_simd->vector4D_Distance(&in4d1,&in4d2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );


		result = resultado1 == resultado2 ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector4D_Distance  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );

}

static void TestVector4DAlignedSum()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	ALIGNTO16 CVec4D tst,dst;
	const char *result;
	in4d1Aligned.toZero();
	in4d2Aligned.toZero();
	reset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector4D Aligned Sum()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector4D_AlignedSum(&in4d1Aligned, &in4d2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = in4d1Aligned;

		CProcClock::printClocks( va( "generic->vector4D_AlignedSum "), 1, bestClocksGeneric );
		reset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector4D_AlignedSum(&in4d1Aligned, &in4d2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = in4d1Aligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector4D_AlignedSum  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestVector4DAlignedSumOf()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	ALIGNTO16 CVec4D tst,dst;
	const char *result;
	in4d1Aligned.toZero();
	in4d2Aligned.toZero();
	reset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector4D Aligned SumOf()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector4D_AlignedSumOf(&out4dAligned, &in4d1Aligned, &in4d2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = out4dAligned;

		CProcClock::printClocks( va( "generic->vector4D_AlignedSumOf "), 1, bestClocksGeneric );
		reset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector4D_AlignedSumOf(&out4dAligned, &in4d1Aligned, &in4d2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = out4dAligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector4D_AlignedSumOf  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestVector4DAlignedDiff()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	ALIGNTO16 CVec4D tst,dst;
	const char *result;
	in4d1Aligned.toZero();
	in4d2Aligned.toZero();
	reset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector4D Aligned Diff()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector4D_AlignedDiff(&in4d1Aligned, &in4d2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = in4d1Aligned;

		CProcClock::printClocks( va( "generic->vector4D_AlignedDiff "), 1, bestClocksGeneric );
		reset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector4D_AlignedDiff(&in4d1Aligned, &in4d2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = in4d1Aligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector4D_AlignedDiff  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestVector4DAlignedDiffOf()
{

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	ALIGNTO16 CVec4D tst,dst;
	const char *result;
	in4d1Aligned.toZero();
	in4d2Aligned.toZero();
	reset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector4D Aligned DiffOf()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector4D_AlignedDiffOf(&out4dAligned, &in4d1Aligned, &in4d2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = out4dAligned;

		CProcClock::printClocks( va( "generic->vector4D_AlignedDiffOf "), 1, bestClocksGeneric );
		reset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector4D_AlignedDiffOf(&out4dAligned, &in4d1Aligned, &in4d2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = out4dAligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector4D_AlignedDiffOf  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}
static void TestVector4DAlignedScale()
{

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	ALIGNTO16 CVec4D tst,dst;
	const char *result;
	in4d1Aligned.toZero();
	in4d2Aligned.toZero();
	reset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector4D Aligned scale()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector4D_AlignedScale(&in4d1Aligned, 100.0f);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = in4d1Aligned;

		CProcClock::printClocks( va( "generic->vector4D_AlignedScale "), 1, bestClocksGeneric );
		reset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector4D_AlignedScale(&in4d1Aligned, 100.0f);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = in4d1Aligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector4D_AlignedScale  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestVector4DAlignedScaleOf()
{

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	ALIGNTO16 CVec4D tst,dst;
	const char *result;
	in4d1Aligned.toZero();
	in4d2Aligned.toZero();
	reset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector4D Aligned ScaleOf()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector4D_AlignedScaleOf(&out4dAligned, &in4d1Aligned, 100.0f);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = out4dAligned;

		CProcClock::printClocks( va( "generic->vector4D_AlignedScaleOf "), 1, bestClocksGeneric );
		reset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector4D_AlignedScaleOf(&out4dAligned, &in4d1Aligned, 100.0f);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = out4dAligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector4D_AlignedScaleOf  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}



static void TestVector4DAlignedDot()
{
	ALIGNTO16 CVec4D in4d1Aligned( 1.0f, 2.0f, 3.0f, 4.0f );
	ALIGNTO16 CVec4D in4d2Aligned( 5.0f, 6.0f, 7.0f, 8.0f );

	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGNTO16 CVec4D tst,dst;
	const char *result;
	in4d1Aligned.toZero();
	in4d2Aligned.toZero();
	reset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector4D Aligned dot()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic-> vector4D_AlignedDot(&in4d1Aligned, &in4d2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = in4d1Aligned;

		CProcClock::printClocks( va( "generic->vector4D_AlignedDot "), 1, bestClocksGeneric );
		reset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector4D_AlignedDot(&in4d1Aligned, &in4d2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = in4d1Aligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector4D_AlignedDot  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );


	/* Should be 5 + 12 + 21 + 32 = 70 */


}


static void TestCVec4DAlignedLengthSq()
{

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	const char *result;
	float resultado1,resultado2;
	in4d1Aligned.toZero();
	in4d2Aligned.toZero();
	reset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector4D Aligned getLengthSqr()..."<<endl;


			StartRecordTimeLocal( start );
			resultado1 = p_generic->vector4D_AlignedLengthSq(&in4d2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		CProcClock::printClocks( va( "generic->vector4D_AlignedLengthSq "), 1, bestClocksGeneric );
		reset4D();
		bestClocksSIMD = 0;


			StartRecordTimeLocal( start );
			resultado2 = p_simd->vector4D_AlignedLengthSq(&in4d2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );


		result = resultado1 == resultado2 ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector4D_AlignedLengthSq  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );

}

static void TestCVec4DAlignedLength()
{
	reset4D();

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	const char *result;
	float resultado1,resultado2;
	in4d1Aligned.toZero();
	in4d2Aligned.toZero();
	reset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector4D Aligned getLenght()..."<<endl;


			StartRecordTimeLocal( start );
			resultado1 = p_generic->vector4D_AlignedLength(&in4d2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );


		CProcClock::printClocks( va( "generic->vector4D_AlignedLength "), 1, bestClocksGeneric );
		reset4D();
		bestClocksSIMD = 0;


			StartRecordTimeLocal( start );
			resultado2=p_simd->vector4D_AlignedLength(&in4d2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );


		result = resultado1 == resultado2 ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector4D_AlignedLength  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );

}


static void TestCVec4DAlignedNormalize()
{

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	ALIGNTO16 CVec4D tst,dst;
	const char *result;
	in4d1Aligned.toZero();
	in4d2Aligned.toZero();
	reset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector4D  Aligned toNormal()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector4D_AlignedNormalize(&in4d1Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = in4d1Aligned;

		CProcClock::printClocks( va( "generic->vector4D_AlignedNormalize "), 1, bestClocksGeneric );
		reset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector4D_AlignedNormalize(&in4d1Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = in4d1Aligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector4D_AlignedNormalize  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



	;
}

static void TestCVec4DAlignedNormalizeOf()
{

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	ALIGNTO16 CVec4D tst,dst;
	const char *result;
	in4d1Aligned.toZero();
	in4d2Aligned.toZero();
	reset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector4D Aligned NormalizeOf()..."<<endl;


			StartRecordTimeLocal( start );
			p_generic->vector4D_AlignedNormalizeOf(&out4dAligned, &in4d1Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

		tst = out4dAligned;

		CProcClock::printClocks( va( "generic->vector4D_AlignedNormalizeOf "), 1, bestClocksGeneric );
		reset4D();
		bestClocksSIMD = 0;

		    dst.toZero();
			StartRecordTimeLocal( start );
			p_simd->vector4D_AlignedNormalizeOf(&out4dAligned, &in4d1Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
		dst = out4dAligned;

		result = dst.compare( tst, MATX_SIMD_EPSILON ) ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector4D_AlignedNormalizeOf  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );



}

static void TestCVec4DAlignedDistance()
{
	reset4D();

	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	const char *result;
	float resultado1,resultado2;
	in4d1Aligned.toZero();
	in4d2Aligned.toZero();
	reset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Vector4D Aligned distance()..."<<endl;


			StartRecordTimeLocal( start );
			resultado1 = p_generic->vector4D_AlignedDistance(&in4d1Aligned,&in4d2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );


		CProcClock::printClocks( va( "generic->vector4D_AlignedDistance "), 1, bestClocksGeneric );
		reset4D();
		bestClocksSIMD = 0;


			StartRecordTimeLocal( start );
			resultado2=p_simd->vector4D_AlignedDistance(&in4d1Aligned,&in4d2Aligned);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );


		result = resultado1 == resultado2 ? "ok" : S_COLOR_RED"X";
		CProcClock::printClocks( va( "   simd->vector4D_AlignedDistance  %s",  result ), 1, bestClocksSIMD, bestClocksGeneric );

}
/*
============
CSIMD::Test4D_f
============
*/
void CSIMD::Test4D_f( const Util::CCMDLineArgs &args ){

	#ifdef _WIN32
	SetThreadPriority( GetCurrentThread(), THREAD_PRIORITY_TIME_CRITICAL );
#endif /* _WIN32 */

	p_simd = processor;
	p_generic = generic;

	if ( CMyString::length( args.Argv( 1 ) ) != 0 ) {
		cpuid_t cpuid = System::getProcessorId();
		CMyString argString = args.Args();

		argString.replace( " ", "" );

		if ( CMyString::compareInsen( argString, "MMX" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) ) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX\n" <<endl;
				return;
			}
			p_simd = new CSIMD_MMX;
			#if !defined (__GNUC__)
		} else if ( CMyString::compareInsen( argString, "3DNow" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) || !( cpuid & CPUID_3DNOW ) ) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX & 3DNow\n"  <<endl;
				return;
			}
			p_simd = new CSIMD_3DNow;
		#endif
		} else if ( CMyString::compareInsen( argString, "SSE" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) || !( cpuid & CPUID_SSE ) ) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX & SSE\n" <<endl;
				return;
			}
			p_simd = new CSIMD_SSE;
		} else if ( CMyString::compareInsen( argString, "SSE2" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) || !( cpuid & CPUID_SSE ) || !( cpuid & CPUID_SSE2 ) ) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX & SSE & SSE2\n"  <<endl;
				return;
			}
			p_simd = new CSIMD_SSE2;
		} else if ( CMyString::compareInsen( argString, "SSE3" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) || !( cpuid & CPUID_SSE ) || !( cpuid & CPUID_SSE2 ) || !( cpuid & CPUID_SSE3 ) ) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX & SSE & SSE2 & SSE3\n"  <<endl;
				return;
			}
			p_simd = new CSIMD_SSE3();
		}else if ( CMyString::compareInsen( argString, "SSE41" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) || !( cpuid & CPUID_SSE ) || !( cpuid & CPUID_SSE2 ) || !( cpuid & CPUID_SSE3 )|| !( cpuid & CPUID_SSE41 ) ) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX & SSE & SSE2 & SSE3 & SSE41\n"  <<endl;
				return;
			}
			p_simd = new CSIMD_SSE3(CPUID_SSE41);//CCPUID_SSE41();
		}else if ( CMyString::compareInsen( argString, "SSE42" ) == 0 ) {
			if ( !( cpuid & CPUID_MMX ) || !( cpuid & CPUID_SSE ) || !( cpuid & CPUID_SSE2 ) || !( cpuid & CPUID_SSE3 )|| !( cpuid & CPUID_SSE41 ) || !( cpuid & CPUID_SSE42 )) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support MMX & SSE & SSE2 & SSE3 & SSE41 & SSE41\n"  <<endl;
				return;
			}
			p_simd = new CSIMD_SSE3(CPUID_SSE42);//CSIMD_SSE42();
		} /*else if ( CMyString::compareInsen( argString, "AltiVec" ) == 0 ) {
			if ( !( cpuid & CPUID_ALTIVEC ) ) {
				Debug::debug(Debug::math,__FUNCTION__) << "CPU does not support AltiVec\n"  <<endl;
				return;
			}
			p_simd = new CSIMD_AltiVec();
		} */else {
			Debug::debug(Debug::math,__FUNCTION__) << "invalid argument, use: MMX, 3DNow, SSE, SSE2, SSE3, AltiVec\n"  <<endl;
			return;
		}
	}

	//idLib::common->SetRefreshOnPrint( true );

	Debug::debug(Debug::math,__FUNCTION__) << "using %s for SIMD processing: "<< p_simd->getName() <<endl;

	CProcClock::getBaseClocks();
	//testes=============================================

	TestVector4DSum();
	TestVector4DSumOf();
	TestVector4DDiff();
	TestVector4DDiffOf();
	TestVector4DScale();
	TestVector4DScaleOf();
	TestVector4DDot();
	TestCVec4DLengthSq();
	TestCVec4DLength();
	TestCVec4DNormalize();
	TestCVec4DNormalizeOf();
	TestCVec4DDistance();
	TestVector4DAlignedSum();
	TestVector4DAlignedSumOf();
	TestVector4DAlignedDiff();
	TestVector4DAlignedDiffOf();
	TestVector4DAlignedScale();
	TestVector4DAlignedScaleOf();
	TestVector4DAlignedDot();
	TestCVec4DAlignedLengthSq();
	TestCVec4DAlignedLength();
	TestCVec4DAlignedNormalize();
	TestCVec4DAlignedNormalizeOf();
	TestCVec4DAlignedDistance();
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


//===============TRIGONOMETRY=====================================

void TestSin() {
	int i;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( float fdst0[361] );
	ALIGN16( float fdst1[361] );
	ALIGN16( float fsrc0[361] );
	ALIGN16( float fsrc1[361] );

	float fdst00[90][4] ;
	float fdst11[90][4] ;
	float fsrc00[90][4] ;
	float fsrc11[90][4] ;

	const char *result;


	for ( i = 0; i < 360; i++ ) {
		fsrc0[i] = DEG2RAD(i);
		fsrc1[i] = DEG2RAD(i);
	}
	int contador =0;
	for ( i = 0; i < 90; i++ ) {

		fsrc00[i][0] = DEG2RAD(contador );
		fsrc11[i][0] = DEG2RAD(contador );
		contador ++;
        fsrc00[i][1] = DEG2RAD(contador );
		fsrc11[i][1] = DEG2RAD(contador );
		contador  ++;
		fsrc00[i][2] = DEG2RAD(contador );
		fsrc11[i][2] = DEG2RAD(contador );
		contador  ++;
		fsrc00[i][3] = DEG2RAD(contador );
		fsrc11[i][3] = DEG2RAD(contador );
		contador  ++;
	}

	Debug::debug(Debug::math,__FUNCTION__) <<"====================================" <<endl;

	bestClocksGeneric = 0;
	for ( i = 0; i < 360; i++ ) {
		StartRecordTimeLocal( start );
		fdst0[i]=p_generic->sin( fsrc0[i] );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->sin( float )", 360, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < 360; i++ ) {
		StartRecordTimeLocal( start );
		fdst1[i]=p_simd->sin( fsrc1[i] );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < 360; i++ ) {
		if ( CMath::fabs( fdst0[i] - fdst1[i] ) > 1e-5f ) {
			break;
		}
	}
	result = ( i >= 360 ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->sin( float ) %s", result ), 360, bestClocksSIMD, bestClocksGeneric );
//=====================================================sinZeroHalfPI

	bestClocksGeneric = 0;
	for ( i = 0; i < 90; i++ ) {
		StartRecordTimeLocal( start );
		fdst0[i]=p_generic->sinZeroHalfPI( fsrc0[i] );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->sinZeroHalfPI( float )", 90, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < 90; i++ ) {
		StartRecordTimeLocal( start );
		fdst1[i]=p_simd->sinZeroHalfPI( fsrc1[i] );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < 90; i++ ) {
		if ( CMath::fabs( fdst0[i] - fdst1[i] ) > 1e-5f ) {
			break;
		}
	}
	result = ( i >= 90 ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->sinZeroHalfPI( float ) %s", result ), 90, bestClocksSIMD, bestClocksGeneric );




//==============================================
	bestClocksGeneric = 0;
	for ( i = 0; i < 90; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->sin4( fsrc00[i],fdst00[i] );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->sin4( float[]  )", 90, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < 90; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->sin4( fsrc11[i],fdst11[i] );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < 90; i++ ) {
		if ( CMath::fabs( fdst00[i][0] - fdst11[i][0] ) > 1e-5f ) {
			break;
		}
		if ( CMath::fabs( fdst00[i][1] - fdst11[i][1] ) > 1e-5f ) {
			break;
		}
		if ( CMath::fabs( fdst00[i][2] - fdst11[i][2] ) > 1e-5f ) {
			break;
		}
		if ( CMath::fabs( fdst00[i][3] - fdst11[i][3] ) > 1e-5f ) {
			break;
		}
	}
	result = ( i >= 90) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->sin4( float[] ) %s", result ), 90, bestClocksSIMD, bestClocksGeneric );
//=========================================================================

bestClocksGeneric = 0;
	for ( i = 0; i < 22; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->sin4ZeroHalfPI( fsrc00[i],fdst00[i] );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->sin4ZeroHalfPI( float[]  )", 22, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < 22; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->sin4ZeroHalfPI( fsrc11[i],fdst11[i] );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < 22; i++ ) {
		if ( CMath::fabs( fdst00[i][0] - fdst11[i][0] ) > 1e-5f ) {
			break;
		}
		if ( CMath::fabs( fdst00[i][1] - fdst11[i][1] ) > 1e-5f ) {
			break;
		}
		if ( CMath::fabs( fdst00[i][2] - fdst11[i][2] ) > 1e-5f ) {
			break;
		}
		if ( CMath::fabs( fdst00[i][3] - fdst11[i][3] ) > 1e-5f ) {
			break;
		}
	}
	result = ( i >= 22) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->sin4ZeroHalfPI( float[] ) %s", result ), 22, bestClocksSIMD, bestClocksGeneric );


}

void TestCos() {
	int i;
	TIME_TYPE start, end, bestClocksGeneric, bestClocksSIMD;
	ALIGN16( float fdst0[361] );
	ALIGN16( float fdst1[361] );
	ALIGN16( float fsrc0[361] );
	ALIGN16( float fsrc1[361] );

	float fdst00[90][4] ;
	float fdst11[90][4] ;
	float fsrc00[90][4] ;
	float fsrc11[90][4] ;

	const char *result;


	for ( i = 0; i < 360; i++ ) {
		fsrc0[i] = DEG2RAD(i);
		fsrc1[i] = DEG2RAD(i);
	}
	int contador =0;
	for ( i = 0; i < 90; i++ ) {

		fsrc00[i][0] = DEG2RAD(contador );
		fsrc11[i][0] = DEG2RAD(contador );
		contador ++;
        fsrc00[i][1] = DEG2RAD(contador );
		fsrc11[i][1] = DEG2RAD(contador );
		contador  ++;
		fsrc00[i][2] = DEG2RAD(contador );
		fsrc11[i][2] = DEG2RAD(contador );
		contador  ++;
		fsrc00[i][3] = DEG2RAD(contador );
		fsrc11[i][3] = DEG2RAD(contador );
		contador  ++;
	}

	Debug::debug(Debug::math,__FUNCTION__) <<"====================================" <<endl;

	bestClocksGeneric = 0;
	for ( i = 0; i < 360; i++ ) {
		StartRecordTimeLocal( start );
		fdst0[i]=p_generic->cos( fsrc0[i] );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->cos( float )", 360, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < 360; i++ ) {
		StartRecordTimeLocal( start );
		fdst1[i]=p_simd->cos( fsrc1[i] );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < 360; i++ ) {
		if ( CMath::fabs( fdst0[i] - fdst1[i] ) > 1e-5f ) {
			break;
		}
	}
	result = ( i >= 360 ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->cos( float ) %s", result ), 360, bestClocksSIMD, bestClocksGeneric );
//=====================================================cosZeroHalfPI

	bestClocksGeneric = 0;
	for ( i = 0; i < 90; i++ ) {
		StartRecordTimeLocal( start );
		fdst0[i]=p_generic->cosZeroHalfPI( fsrc0[i] );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->cosZeroHalfPI( float )", 90, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < 90; i++ ) {
		StartRecordTimeLocal( start );
		fdst1[i]=p_simd->cosZeroHalfPI( fsrc1[i] );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < 90; i++ ) {
		if ( CMath::fabs( fdst0[i] - fdst1[i] ) > 1e-5f ) {
			break;
		}
	}
	result = ( i >= 90 ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->cosZeroHalfPI( float ) %s", result ), 90, bestClocksSIMD, bestClocksGeneric );




//==============================================
	bestClocksGeneric = 0;
	for ( i = 0; i < 90; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->cos4( fsrc00[i],fdst00[i] );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->cos4( float[]  )", 90, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < 90; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->cos4( fsrc11[i],fdst11[i] );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < 90; i++ ) {
		if ( CMath::fabs( fdst00[i][0] - fdst11[i][0] ) > 1e-5f ) {
			break;
		}
		if ( CMath::fabs( fdst00[i][1] - fdst11[i][1] ) > 1e-5f ) {
			break;
		}
		if ( CMath::fabs( fdst00[i][2] - fdst11[i][2] ) > 1e-5f ) {
			break;
		}
		if ( CMath::fabs( fdst00[i][3] - fdst11[i][3] ) > 1e-5f ) {
			break;
		}
	}
	result = ( i >= 90) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->cos4( float[] ) %s", result ), 90, bestClocksSIMD, bestClocksGeneric );
//=========================================================================

bestClocksGeneric = 0;
	for ( i = 0; i < 22; i++ ) {
		StartRecordTimeLocal( start );
		p_generic->cos4ZeroHalfPI( fsrc00[i],fdst00[i] );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksGeneric );
	}
	CProcClock::printClocks( "generic->cos4ZeroHalfPI( float[]  )", 22, bestClocksGeneric );

	bestClocksSIMD = 0;
	for ( i = 0; i < 22; i++ ) {
		StartRecordTimeLocal( start );
		p_simd->cos4ZeroHalfPI( fsrc11[i],fdst11[i] );
		StopRecordTimeLocal( end );
		GetBest( start, end, bestClocksSIMD );
	}

	for ( i = 0; i < 22; i++ ) {
		if ( CMath::fabs( fdst00[i][0] - fdst11[i][0] ) > 1e-5f ) {
			break;
		}
		if ( CMath::fabs( fdst00[i][1] - fdst11[i][1] ) > 1e-5f ) {
			break;
		}
		if ( CMath::fabs( fdst00[i][2] - fdst11[i][2] ) > 1e-5f ) {
			break;
		}
		if ( CMath::fabs( fdst00[i][3] - fdst11[i][3] ) > 1e-5f ) {
			break;
		}
	}
	result = ( i >= 22) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->cos4ZeroHalfPI( float[] ) %s", result ), 22, bestClocksSIMD, bestClocksGeneric );


}

/*
============
SSE_TestTrigonometry
============
*/
void CSIMD::SSE_TestTrigonometry() {

#ifdef _WIN32
	SetThreadPriority( GetCurrentThread(), THREAD_PRIORITY_TIME_CRITICAL );
#endif /* _WIN32 */

	p_simd = processor;
	p_generic = generic;


	TestSin();
	TestCos();
//old ones
#if 0
	int i;
	float a, s1, s2, c1, c2;

	for ( i = 0; i < 100; i++ ) {
		a = i * CMath::HALF_PI / 100.0f;

		s1 = sin( a );
		s2 = p_simd->sinZeroHalfPI( a );

		if ( fabs( s1 - s2 ) > 1e-7f ) {
			SMF_ASSERT( 0 );
		}

		c1 = CMath::cos( a );
		c2 = p_simd->cosZeroHalfPI( a );

		if ( fabs( c1 - c2 ) > 1e-7f ) {
			SMF_ASSERT( 0 );
		}
	}

	for ( i = -200; i < 200; i++ ) {
		a = i * CMath::TWO_PI / 100.0f;

		s1 = sin( a );
		s2 = p_simd->sin( a );

		if ( fabs( s1 - s2 ) > 1e-6f ) {
			SMF_ASSERT( 0 );
		}

		c1 = CMath::cos( a );
		c2 = p_simd->cos( a );

		if ( fabs( c1 - c2 ) > 1e-6f ) {
			SMF_ASSERT( 0 );
		}

	//	SSE_SinCos( a, s2, c2 );
	//	if ( fabs( s1 - s2 ) > 1e-6f || fabs( c1 - c2 ) > 1e-6f ) {
	//		SMF_ASSERT( 0 );
	//	}
	}
#endif
    if ( p_simd != processor ) {
		delete p_simd;
	}
	p_simd = NULL;
	p_generic = NULL;

#ifdef _WIN32
	SetThreadPriority( GetCurrentThread(), THREAD_PRIORITY_NORMAL );
#endif /* _WIN32 */
}


//==============END TRIGINOMETRY=============================
} //END MATH
} //end SMF
