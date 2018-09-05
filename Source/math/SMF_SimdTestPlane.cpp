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
#include "geometry/SMF_Plane.h"
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




//===================CQuaternion====================================



CPlane planout, planin1, planin2;
ALIGNTO16  CPlane planoutAligned;
ALIGNTO16  CPlane planin1Aligned;
ALIGNTO16  CPlane planin2Aligned;
static void matReset4D()
{
	planout.toZero();
	planoutAligned.toZero();


	planin1.a=1.0f;
	planin1.b=2.0f;
	planin1.c=3.0f;
	planin1.d=4.0f;

	planin2.a=40.0f;
	planin2.b=30.0f;
	planin2.c=2.0f;
	planin2.d=10.0f;

	planin1Aligned = planin1;
	planin2Aligned = planin2;

}



static void Test_plane_FromPoints()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CPlane tst[NUMTESTS],dst[NUMTESTS];
	CVec3D pa[NUMTESTS],pb[NUMTESTS],pc[NUMTESTS];
	const char *result;
	matReset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing  plane_FromPoints..."<<endl;

	CRandom srnd( RANDOM_SEED );

	for ( int i = 0; i < NUMTESTS; i++ ) {
		dst[i].toZero();
		tst[i].toZero();
		pa[i]= CVec3D(srnd.randomFloat(),srnd.randomFloat(),srnd.randomFloat());
		pb[i]= CVec3D(srnd.randomFloat(),srnd.randomFloat(),srnd.randomFloat());
		pc[i]= CVec3D(srnd.randomFloat(),srnd.randomFloat(),srnd.randomFloat());
	}
	        int i=0;
			for ( i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			p_generic-> plane_FromPoints(&planout,&pa[i],&pb[i],&pc[i]);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
			tst[i] = planout;
			}


		CProcClock::printClocks( va( "generic-> plane_FromPoints "), NUMTESTS, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;



			for ( i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			p_simd-> plane_FromPoints(&planout,&pa[i],&pb[i],&pc[i]);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
			dst[i] = planout;
			}
		for ( i = 0; i < NUMTESTS; i++ ) {
			if ( !dst[i].compare( tst[i], MATX_SIMD_EPSILON ) ) {
				break;
			}
		}
	result = ( i >= NUMTESTS ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd-> plane_FromPoints  %s",  result ), NUMTESTS, bestClocksSIMD, bestClocksGeneric );
}


static void Test_plane_DistToPoint()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	float tst[NUMTESTS],dst[NUMTESTS];
	CVec3D pa[NUMTESTS];
	const char *result;
	matReset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing  plane_DistToPoint..."<<endl;

	CRandom srnd( RANDOM_SEED );

	for ( int i = 0; i < NUMTESTS; i++ ) {
		dst[i]=0.0f;
		tst[i]=0.0f;
		pa[i]= CVec3D(srnd.randomFloat(),srnd.randomFloat(),srnd.randomFloat());

	}
	        int i=0;
			for ( i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			tst[i]=p_generic-> plane_DistToPoint(&planin1,&pa[i]);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

			}


		CProcClock::printClocks( va( "generic-> plane_DistToPoint "), NUMTESTS, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;



			for ( i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			dst[i] =p_simd-> plane_DistToPoint(&planin1,&pa[i]);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );

			}
		for ( i = 0; i < NUMTESTS; i++ ) {
			if ( !dst[i] == tst[i] ) {
				break;
			}
		}
	result = ( i >= NUMTESTS ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd-> plane_DistToPoint  %s",  result ), NUMTESTS, bestClocksSIMD, bestClocksGeneric );
}



static void Test_plane_Dot()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	float tst[NUMTESTS],dst[NUMTESTS];
	CVec3D pa[NUMTESTS];
	const char *result;
	matReset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Test_plane_Dot..."<<endl;

	CRandom srnd( RANDOM_SEED );

	for ( int i = 0; i < NUMTESTS; i++ ) {
		dst[i]=0.0f;
		tst[i]=0.0f;
		pa[i]= CVec3D(srnd.randomFloat(),srnd.randomFloat(),srnd.randomFloat());

	}
	        int i=0;
			for ( i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			tst[i]=p_generic-> plane_Dot(&planin1,&pa[i]);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

			}


		CProcClock::printClocks( va( "generic->Test_plane_Dot "), NUMTESTS, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;



			for ( i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			dst[i] =p_simd-> plane_Dot(&planin1,&pa[i]);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );

			}
		for ( i = 0; i < NUMTESTS; i++ ) {
			if ( !dst[i] == tst[i] ) {
				break;
			}
		}
	result = ( i >= NUMTESTS ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->Test_plane_Dot  %s",  result ), NUMTESTS, bestClocksSIMD, bestClocksGeneric );
}


static void Test_plane_Dot4()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	float tst[NUMTESTS],dst[NUMTESTS];
	CVec4D pa[NUMTESTS];
	const char *result;
	matReset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing Test_plane_Dot4..."<<endl;

	CRandom srnd( RANDOM_SEED );

	for ( int i = 0; i < NUMTESTS; i++ ) {
		dst[i]=0.0f;
		tst[i]=0.0f;
		pa[i]= CVec4D(srnd.randomFloat(),srnd.randomFloat(),srnd.randomFloat(),srnd.randomFloat());

	}
	        int i=0;
			for ( i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			tst[i]=p_generic-> plane_Dot4(&planin1,&pa[i]);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

			}


		CProcClock::printClocks( va( "generic->Test_plane_Dot4 "), NUMTESTS, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;



			for ( i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			dst[i] =p_simd-> plane_Dot4(&planin1,&pa[i]);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );

			}
		for ( i = 0; i < NUMTESTS; i++ ) {
			if ( !dst[i] == tst[i] ) {
				break;
			}
		}
	result = ( i >= NUMTESTS ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd->Test_plane_Dot4  %s",  result ), NUMTESTS, bestClocksSIMD, bestClocksGeneric );
}


static void Test_plane_DotNormal()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	float tst[NUMTESTS],dst[NUMTESTS];
	CPlane pa[NUMTESTS];
	CVec3D pnormal[NUMTESTS];
	const char *result;
	matReset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing  plane_DotNormal..."<<endl;

	CRandom srnd( RANDOM_SEED );

	for ( int i = 0; i < NUMTESTS; i++ ) {
		dst[i]=0.0f;
		tst[i]=0.0f;
		pa[i]= CPlane(srnd.randomFloat(),srnd.randomFloat(),srnd.randomFloat(),srnd.randomFloat());
		pnormal[i]=pa[i].getNormal();
	}
	        int i=0;
			for ( i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			tst[i]=p_generic-> plane_DotNormal(&planin1,&pnormal[i]);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

			}


		CProcClock::printClocks( va( "generic-> plane_DotNormal "), NUMTESTS, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;



			for ( i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			dst[i] =p_simd-> plane_DotNormal(&planin1,&pnormal[i]);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );

			}
		for ( i = 0; i < NUMTESTS; i++ ) {
			if ( !dst[i] == tst[i] ) {
				break;
			}
		}
	result = ( i >= NUMTESTS ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd-> plane_DotNormal  %s",  result ), NUMTESTS, bestClocksSIMD, bestClocksGeneric );
}


static void Test_plane_DotPlane()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	float tst[NUMTESTS],dst[NUMTESTS];
	const char *result;
	matReset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing  plane_DotPlane..."<<endl;

	CRandom srnd( RANDOM_SEED );

	for ( int i = 0; i < NUMTESTS; i++ ) {
		dst[i]=0.0f;
		tst[i]=0.0f;

	}
	        int i=0;
			for ( i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			tst[i]=p_generic-> plane_DotPlane(&planin1,&planin2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );

			}


		CProcClock::printClocks( va( "generic-> plane_DotPlane "), NUMTESTS, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;



			for ( i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			dst[i] =p_simd-> plane_DotPlane(&planin1,&planin2);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );

			}
		for ( i = 0; i < NUMTESTS; i++ ) {
			if ( !dst[i] == tst[i] ) {
				break;
			}
		}
	result = ( i >= NUMTESTS ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd-> plane_DotPlane  %s",  result ), NUMTESTS, bestClocksSIMD, bestClocksGeneric );
}


static void Test_plane_Normalize()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CPlane tst[NUMTESTS],dst[NUMTESTS];
	const char *result;
	matReset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing  plane_Normalize..."<<endl;

	CRandom srnd( RANDOM_SEED );

	for ( int i = 0; i < NUMTESTS; i++ ) {
		dst[i].toZero();
		tst[i].toZero();
	}
	        int i=0;
			for ( i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			p_generic-> plane_Normalize(&planin1);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
			tst[i] = planin1;
			}


		CProcClock::printClocks( va( "generic-> plane_Normalize "), NUMTESTS, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;



			for ( i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			p_simd-> plane_Normalize(&planin1);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
			dst[i] = planin1;
			}
		for ( i = 0; i < NUMTESTS; i++ ) {
			if ( !dst[i].compare( tst[i], MATX_SIMD_EPSILON ) ) {
				break;
			}
		}
	result = ( i >= NUMTESTS ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd-> plane_Normalize  %s",  result ), NUMTESTS, bestClocksSIMD, bestClocksGeneric );
}

static void Test_plane_NormalizeOf()
{
	TIME_TYPE start=0, end=0, bestClocksGeneric=0, bestClocksSIMD=0;
	CPlane tst[NUMTESTS],dst[NUMTESTS];
	const char *result;
	matReset4D();
	Debug::debug(Debug::math,__FUNCTION__) << "Testing  plane_NormalizeOf..."<<endl;

	CRandom srnd( RANDOM_SEED );

	for ( int i = 0; i < NUMTESTS; i++ ) {
		dst[i].toZero();
		tst[i].toZero();
	}
	        int i=0;
			for ( i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			p_generic-> plane_NormalizeOf(&planout,&planin1);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksGeneric );
			tst[i] = planout;
			}


		CProcClock::printClocks( va( "generic-> plane_NormalizeOf "), NUMTESTS, bestClocksGeneric );
		matReset4D();
		bestClocksSIMD = 0;



			for ( i = 0; i < NUMTESTS; i++ ) {
			StartRecordTimeLocal( start );
			p_simd-> plane_NormalizeOf(&planout,&planin1);
			StopRecordTimeLocal( end );
			GetBest( start, end, bestClocksSIMD );
			dst[i] = planout;
			}
		for ( i = 0; i < NUMTESTS; i++ ) {
			if ( !dst[i].compare( tst[i], MATX_SIMD_EPSILON ) ) {
				break;
			}
		}
	result = ( i >= NUMTESTS ) ? "ok" : S_COLOR_RED"X";
	CProcClock::printClocks( va( "   simd-> plane_NormalizeOf  %s",  result ), NUMTESTS, bestClocksSIMD, bestClocksGeneric );
}

void CSIMD::TestPlane4D_f(){

#ifdef _WIN32
	SetThreadPriority( GetCurrentThread(), THREAD_PRIORITY_TIME_CRITICAL );
#endif /* _WIN32 */

	p_simd = processor;
	p_generic = generic;



	//idLib::common->SetRefreshOnPrint( true );

	Debug::debug(Debug::math,__FUNCTION__) << "using %s for SIMD processing: "<< p_simd->getName() <<endl;

	CProcClock::getBaseClocks();
	//testes=============================================

Test_plane_FromPoints();
Test_plane_DistToPoint();
Test_plane_Dot();
Test_plane_Dot4();
Test_plane_DotNormal();
Test_plane_DotPlane();
Test_plane_Normalize();
Test_plane_NormalizeOf();

#if 0

#endif

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
