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


#ifndef _SMF__SYSTEM_MINGW__
#define _SMF__SYSTEM_MINGW__


/* system utilities */
#include <stdint.h>

#include <string>
#include "../SMF_Config.h"
#include "../util/SMF_UtilStructs.h"
#include "../math/SMF_Math.h"
namespace SMF {
namespace System {

/* Here is the dirty part. Set up your 387 through the control word
 * (cw) register.
 *
 *     15-13    12  11-10  9-8     7-6     5    4    3    2    1    0
 * | reserved | IC | RC  | PC | reserved | PM | UM | OM | ZM | DM | IM
 *
 * IM: Invalid operation mask
 * DM: Denormalized operand mask
 * ZM: Zero-divide mask
 * OM: Overflow mask
 * UM: Underflow mask
 * PM: Precision (inexact result) mask
 *
 * Mask bit is 1 means no interrupt.
 *
 * PC: Precision control
 * 11 - round to extended precision
 * 10 - round to double precision
 * 00 - round to single precision
 *
 * RC: Rounding control
 * 00 - rounding to nearest
 * 01 - rounding down (toward - infinity)
 * 10 - rounding up (toward + infinity)
 * 11 - rounding toward zero
 *
 * IC: Infinity control
 * That is for 8087 and 80287 only.
 *
 * The hardware default is 0x037f which we use.
 */



/** \brief FPU rounding control 
\warning To use directly with fpu_SetRounding. Not to use directly on the fpu control word, if you need to set up directly the control word use the defines
\see fpu_SetRounding
**/
typedef enum {
	FPU_ROUNDING_TO_NEAREST				= 0x0, /* RECOMMENDED */
	FPU_ROUNDING_DOWN					= 1,
	FPU_ROUNDING_UP						= 2,
	FPU_ROUNDING_TO_ZERO				= 3
} fpuRounding_t;
/* rounding control */
/*words to use directly on the FPU control register rounding control 
*/
#define _FPU_RC_NEAREST 0x0    /* RECOMMENDED */
#define _FPU_RC_DOWN    0x400
#define _FPU_RC_UP      0x800
#define _FPU_RC_ZERO    0xC00
#define _FPU_RESERVED 0xF0C0  /* Reserved bits in cw */

/** 
\brief defines FPU Control bits
\note words can be used directly on the FPU control word
FPU_EXCEPTION_INVALID_OPERATION   0x01
FPU_EXCEPTION_DENORMALIZED_OPERAND 0x02
FPU_EXCEPTION_DIVIDE_BY_ZERO   0x04
FPU_EXCEPTION_NUMERIC_OVERFLOW   0x08
FPU_EXCEPTION_NUMERIC_UNDERFLOW   0x10
FPU_EXCEPTION_INEXACT_RESULT   0x20
FPU_EXCEPTION_ALL   (FE_INEXACT | FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID)
\see fpu_EnableExceptions
**/
typedef enum {
	FPU_EXCEPTION_INVALID_OPERATION		= 0x01,
	FPU_EXCEPTION_DENORMALIZED_OPERAND	= 0x02,
	FPU_EXCEPTION_DIVIDE_BY_ZERO		= 0x04,
	FPU_EXCEPTION_NUMERIC_OVERFLOW		= 0x08,
	FPU_EXCEPTION_NUMERIC_UNDERFLOW		= 0x10,
	FPU_EXCEPTION_INEXACT_RESULT		= 0x20,
	FPU_EXCEPTION_ALL					= (FPU_EXCEPTION_INVALID_OPERATION || FPU_EXCEPTION_DENORMALIZED_OPERAND || FPU_EXCEPTION_DIVIDE_BY_ZERO || FPU_EXCEPTION_NUMERIC_OVERFLOW || FPU_EXCEPTION_NUMERIC_UNDERFLOW || FPU_EXCEPTION_INEXACT_RESULT)
} fpuExceptions_t;

/** 
\brief precision control to use with the fpuSetPrecision function
\warning do not use directly on the fpu control word. If you need to setup the fpu control word use the defines (_FPU_EXTENDED, _FPU_DOUBLE, _FPU_SINGLE)
\see fpuSetPrecision
**/
typedef enum {
	FPU_PRECISION_SINGLE				= 0,
	FPU_PRECISION_DOUBLE				= 2,
	FPU_PRECISION_DOUBLE_EXTENDED		= 3
} fpuPrecision_t;

/* precision control words to use directly on the FPU control register*/
#define _FPU_EXTENDED 0x300	/* libm requires double extended precision.  */
#define _FPU_DOUBLE   0x200
#define _FPU_SINGLE   0x0

/* The fdlibm code requires strict IEEE double precision arithmetic,
   and no interrupts for exceptions, rounding to nearest.  */

#define _FPU_DEFAULT  0x037f

/* IEEE:  same as above.  */
#define _FPU_IEEE     0x037f
/// Statistics about the memory use
typedef struct sysMemoryStats_s {
	int memoryLoad;
	int totalPhysical;
	int availPhysical;
	int totalPageFile;
	int availPageFile;
	int totalVirtual;
	int availVirtual;
	int availExtendedVirtual;
} sysMemoryStats_t;

typedef unsigned long address_t;

//template<class type> class CList;		// for Sys_ListFiles
//=========Clock===========================================================

/**
 * \class CTime
 *
 * \ingroup SMF_System
 *
 * \if pt_br
 * \brief high resolution clock and timestamp control 
 * \elseif us_en
 * \brief clock and timestamp control
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class CTime{
	public:
static uint64_t currentMicroseconds();
static unsigned long memoryUsage();
static uint64_t currentSeconds();


// a decent minimum sleep time to avoid going below the process scheduler speeds
#define	SYS_MINSLEEP	20

// Sys_Milliseconds should only be used for profiling purposes,
// any game related timing information should come from event timestamps
// return the elapsed time in mileseconds from the first running of the function.
static int		getElapsedMilliseconds();

/*
\brief get time stamp information from CPU for accurate performance testing
\note example: how to get correct timestamp information:
\n cpuid ; force all previous instructions to complete
\n rdtsc ; read time stamp counter
\n mov time, eax ; move counter into variable
\n fdiv ; floating-point divide (medindo quantos ciclos demora para fazer a divisão)
\n cpuid ; wait for FDIV to complete before RDTSC
\n rdtsc ; read time stamp counter
\n sub eax, time ; find the difference
\note Using the "rdtsc"
\n This method use rdtsc.
\n rdtsc is a time stamp counter that returns  the number of clock ticks from the time the system was last reset.
\n The rdtsc instruction returns the time stamp counter in EDX:EAX.

\n This instruction is useful when you want to measure the performance of a certain code or application in clock ticks,
\n or compare the performance of two programs, which are too small to be counted in seconds.
*/

static double		getClockTicks();
static double		getClockTicksPerSecond();
/**
\brief record the current time stamp information to the start var
\note save the actual ebx state to ::SMF::MATH::saved_ebx
**/
static void startRecordingTimeStamp(double &start);
/**
\brief record the current time stamp information to the start end
\note restore  the  ebx state saved by startRecordingTimeStamp
**/
static void stopRecordingTimStamp(double &end);
};
//=======================================
// FLOAT Point UNIT
//=======================================


// returns a selection of the CPUID_* flags
cpuid_t		getProcessorId();
const char * getProcessorString();


/**
 * \class CFPU
 *
 * \ingroup SMF_System
 *
 * \if pt_br
 * \brief Set up Floating Point Unit 
 * \elseif us_en
 * \brief Set up Floating Point Unit 
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */

class CFPU{
public:
#if defined(__GNUC__)
//=======================================
// fpu_HasDazSupport
//=======================================
// verifica se existe suporte para DAZ
static void fpu_HasDazSupport(bool &result);
//=======================================
// fpu_getFXSave
//=======================================
// executa fxsave e retorna o resultado no buffer

static void fpu_getFXSave (void *buffer);

#endif
// returns true if the FPU stack is empty
static bool		fpu_StackIsEmpty();

// empties the FPU stack
static void		fpu_ClearStack();

// returns the FPU state as a string
static const char * 	fpu_GetState();

/**
\brief  enables the given FPU exceptions
\note
    \n FPU_EXCEPTION_INVALID_OPERATION		= domain error occurred in an earlier floating-point operation
	\n FPU_EXCEPTION_DENORMALIZED_OPERAND	= 
	\n FPU_EXCEPTION_DIVIDE_BY_ZERO			= pole error occurred in an earlier floating-point operation
	\n FPU_EXCEPTION_NUMERIC_OVERFLOW		= the result of the earlier floating-point operation was too large to be representable
	\n FPU_EXCEPTION_NUMERIC_UNDERFLOW		= the result of the earlier floating-point operation was subnormal with a loss of precision
	\n FPU_EXCEPTION_INEXACT_RESULT			= inexact result: rounding was necessary to store the result of an earlier floating-point operation
	\n FPU_EXCEPTION_ALL					= bitwise OR of all supported floating-point exceptions
**/
static void		fpu_EnableExceptions( fpuExceptions_t exceptions );
/**
\brief  disables the given FPU exceptions
// FE_DIVBYZERO, FE_INEXACT, FE_INVALID, FE_OVERFLOW, FE_UNDERFLOW
**/
static void 		fpu_DisableExceptions( int exceptions );

// sets the FPU precision
/*
The PC field (bits 9 and 8) or Precision Control determines to what precision the FPU rounds results after each arithmetic instruction in one of three ways:

00 = 24 bits (REAL4) = SINGLE PRECISION(0X00)
01 = Not used
10 = 53 bits (REAL8)   =  DOUBLE PRECISION (0x200)
11 = 64 bits (REAL10) (this is the initialized state) = EXTENDED (0x300)

\param precision 0 = single precision
                 1 = not used
				 2 = double precision
				 3 = extended precision
*/
static void		fpu_SetPrecision( fpuPrecision_t precision );

// sets the FPU rounding mode


#if defined(WIN32)  && defined(_MSC_VER)

// x87 fpu
#define getx87cr(x)    __asm fnstcw x;
#define setx87cr(x)    __asm fldcw x;
#define getx87sr(x)    __asm fnstsw x;

// SIMD, gcc with Intel Core 2 Duo uses SSE2(4)
#define getmxcsr(x)   __asm stmxcsr   x;
#define setmxcsr(x)    __asm ldmxcsr   x;

#define saveFX(fx,err)  __asm fxsave   state ;
#define loadFX(fx)	__asm __asm ldmxcsr   x;
#else
// x87 fpu
#define getx87cr(x)    asm ("fnstcw %0" : "=m" (x));
#define setx87cr(x)    asm ("fldcw %0"  : "=m" (x));
#define getx87sr(x)    asm ("fnstsw %0" : "=m" (x));
// SIMD, gcc with Intel Core 2 Duo uses SSE2(4)
#define getmxcsr(fx)    asm volatile("stmxcsr %0" : "=m" (fx));
#define setmxcsr(fx)    asm volatile("ldmxcsr %0" : "=m" (fx));

#endif

/* 
\brief set FPU rounding
\note The RC field (bits 11 and 10) or Rounding Control determines how the FPU will round results in one of four ways:

00 = round to nearest, or to even if equidistant (this is the initialized state)
01 = round down (toward -infinity)
10 = round up (toward +infinity)
11 = truncate (toward 0)
\param rounding 0 NEAREST
				1 DOWN
				2 UP
				3 ZERO
*/
static void	fpu_SetRounding( fpuRounding_t rounding );
//TODO: getRouding
static int16_t	fpu_getRounding( int16_t rounding );
/**\todo fpu_restoreRouding*/
//static void fpu_restoreRouding(fpu_control fpu);
/**
\brief set CPUID_FTZ and enable/disable FTZ
\see http://software.intel.com/en-us/articles/x87-and-sse-floating-point-assists-in-ia-32-flush-to-zero-ftz-and-denormals-are-zero-daz
**/
static void    fpu_enableFTZ();
static void 	fpu_disableFTZ();
// sets Flush-To-zero mode (only available when CPUID_FTZ is set)
static void	fpu_SetFTZ( bool enable );

// set CPUID_DAZ and enable/disable DAZ
static void	fpu_enableDAZ();
static void    fpu_disableDAZ();
// sets Denormals-Are-zero mode (only available when CPUID_DAZ is already set)
static void	fpu_SetDAZ( bool enable );
static int fpu_PrintStateFlags( char *ptr, int ctrl, int stat, int tags, int inof, int inse, int opof, int opse );
}; //end CFPU
//================================================
//  END FPU
//=================================================

// lock and unlock memory
/*Locks the specified region of the process's virtual address space into physical memory,
ensuring that subsequent access to the region will not incur a page fault.
*/
bool		lockMemory( void *ptr, int bytes );
bool		unlockMemory( void *ptr, int bytes );

/**Set FPU Mode of working
\param precision
\param rounding
\param exception_mask
*/
SMF_INLINE static int Sys_ieee_set_FPU_mode (int precision, int rounding, int exception_mask)
{
  unsigned short mode = 0 ;

  switch (precision)
    {
    case SMF_IEEE_SINGLE_PRECISION:
      mode |= _FPU_SINGLE ;
      break ;
    case SMF_IEEE_DOUBLE_PRECISION:
      mode |= _FPU_DOUBLE ;
      break ;
    case SMF_IEEE_EXTENDED_PRECISION:
      mode |= _FPU_EXTENDED ;
      break ;
    default:
      mode |= _FPU_EXTENDED ;
    }

  switch (rounding)
    {
    case SMF_IEEE_ROUND_TO_NEAREST:
      mode |= _FPU_RC_NEAREST ;
      break ;
    case SMF_IEEE_ROUND_DOWN:
      mode |= _FPU_RC_DOWN ;
      break ;
    case SMF_IEEE_ROUND_UP:
      mode |= _FPU_RC_UP ;
      break ;
    case SMF_IEEE_ROUND_TO_ZERO:
      mode |= _FPU_RC_ZERO ;
      break ;
    default:
      mode |= _FPU_RC_NEAREST ;
    }

  if (exception_mask & SMF_IEEE_MASK_INVALID)
    mode |= System::FPU_EXCEPTION_INVALID_OPERATION ;

  if (exception_mask & SMF_IEEE_MASK_DENORMALIZED)
    mode |= System::FPU_EXCEPTION_DENORMALIZED_OPERAND ;

  if (exception_mask & SMF_IEEE_MASK_DIVISION_BY_ZERO)
    mode |= System::FPU_EXCEPTION_DIVIDE_BY_ZERO ;

  if (exception_mask & SMF_IEEE_MASK_OVERFLOW)
    mode |= System::FPU_EXCEPTION_NUMERIC_OVERFLOW ;

  if (exception_mask & SMF_IEEE_MASK_UNDERFLOW)
    mode |= System::FPU_EXCEPTION_NUMERIC_UNDERFLOW ;

  if (exception_mask & SMF_IEEE_TRAP_INEXACT)
    {
      mode &= ~ System::FPU_EXCEPTION_INEXACT_RESULT ;
    }
  else
    {
      mode |=  System::FPU_EXCEPTION_INEXACT_RESULT ;
    }

  setx87cr(mode) ;

  return SMF_OK ;
}


} //end Global
} //end SMF
#endif /* !__SYS_PUBLIC__ */
