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
#ifndef _SMF__CONFIG_H
#define _SMF__CONFIG_H

#include "SMF_Types.h"
#include "util/SMF_Debug.h"


////////////////////////////////////////////////////////////
// Identify the operating system
////////////////////////////////////////////////////////////
/**
 *
 *

Alloca allocates memory from Stack rather then heap which is case in malloc. So, when I return from the routine the memory is freed. So, actually this solves my problem of freeing up of dynamically allocated memory . Freeing of memory allocated through malloc is a major headache and if somehow missed leads to all sorts memory problems.

So, my question is that in spite of the above features still alloca use is discouraged, why?
The answer is right there in the man page (at least on Linux):

    RETURN VALUE The alloca() function returns a pointer to the beginning of the allocated space. If the allocation causes stack overflow, program behaviour is undefined.

Which isn't to say it should never be used. One of the OSS projects I work on uses it extensively, and as long as you're not abusing it (alloca'ing huge values), it's fine. Once you go past the "few hundred bytes" mark, it's time to use malloc and friends, instead. You may still get allocation failures, but at least you'll have some indication of the failure instead of just blowing out the stack.
http://stackoverflow.com/questions/1018853/why-is-alloca-not-considered-good-practice
 *
 *
 * **/
#if defined(_MSC_VER) || defined(__GNUC__)

//TODO  discover if it is 32bit or 64bits
    #include <malloc.h>

    // Windows
    #define SMF_ON_WINDOWS
    #ifndef WIN32_LEAN_AND_MEAN
        #define WIN32_LEAN_AND_MEAN
    #endif
    #ifndef NOMINMAX
        #define NOMINMAX
    #endif
    #pragma warning( disable : 4290 )
    #pragma warning( disable : 4800 )
    #pragma warning( disable : 4996 )
    #pragma warning( disable : 4561 )

    #define BUILD_OS_ID					0
    #define	CPUSTRING					"x86"
    #define CPU_EASYARGS					1
    #if defined(__GNUC__)  // MinGW
        #define	BUILD_STRING					"win-x86-MinGw"
        #define ALIGN16( x )					__attribute__((  aligned(16)  )) x
        #define ALIGN32( x )					__attribute__((  aligned(32)  )) x
		#define ALIGNTO16                       __attribute((aligned(16)))
		#define ALIGNTO32						__attribute((aligned(32)))
        #define _alloca							alloca
        #define _alloca16( x )					((void *)((((int)alloca( (x)+15 )) + 15) & ~15))
        #define _allocafloat16( x )				((void *)((((uintptr_t)alloca( (x)+15 )) + 15) & ~15))
        #define PACKED							__attribute__((packed))
		#define id_attribute(x)					__attribute__(x)
		#define VPCALL							__attribute__((fastcall))
		#define CONST_WIN32
#if not defined(FLT_MAX)
#define FLT_MAX 3.40282347e+38F
#endif // not
    #elif defined(_MSC_VER)   //MSVC
        #define	BUILD_STRING					"win-x86"
        #define ALIGN16( x )					__declspec(align(16)) x
        #define ALIGN32( x )					__declspec(align(32)) x
		#define ALIGNTO16                       __declspec(align(16))
        #define ALIGNTO32						__declspec(align(32))
		#define PACKED
        #define _alloca16( x )					((void *)((((int)_alloca( (x)+15 )) + 15) & ~15))
        #define _allocafloat16( x )				((void *)((((uintptr_t)alloca( (x)+15 )) + 15) & ~15))
		#define id_attribute(x)
		#define VPCALL							__fastcall
// The CONST_WIN32 is a #define which resolves to 'const' on Windows, and null on other
// platforms. This #define is used on Windows to detect accidental programming errors
// occurring from an expression "const CVec3D vec; vec[1] = 5;". Trying to return
// const float from operator[] on GCC gives a warning "type qualifiers ignored on function return type",
// so hence this is only enabled on Visual Studio.

		#define CONST_WIN32 const



#endif
    //definitions valid for mingw or msvc
	#define PATHSEPERATOR_STR				"\\"
    #define PATHSEPERATOR_CHAR				'\\'

//=========MEMORY HEAP DEFINES=========================
/// use libc malloc instead of CHeap class one
#ifndef USE_LIBC_MALLOC
	#define USE_LIBC_MALLOC		0
#endif

#ifndef CRASH_ON_STATIC_ALLOCATION
//	#define CRASH_ON_STATIC_ALLOCATION
#endif

//=========ASSERTS AND FAILS============================
// IF NO FAILS DEFINE IS USED THE SYSTEM WILL ONLY LOG THE FAIL

// If MATH_ENABLE_INSECURE_OPTIMIZATIONS is defined, some input data is assumed to be correct and will
// not be checked against at runtime.
// If this flag is undefined (the default), some input is sanity checked so that user cannot crash the system
// e.g. with out-of-bounds accesses.
#define MATH_ENABLE_INSECURE_OPTIMIZATIONS

#define FAIL_USING_EXCEPTIONS
///FAIL_USING_ASSERT WILL CAUSE ASSERT to be raised on fails
#ifdef FAIL_USING_ASSERT
#define SMF_ASSERT(x) SMF_ASSERT(x)
/// FAIL_USING_SILENT WILL CAUSE NOTHING TO BE DONE ON FAILS
#elif defined(FAIL_USING_SILENT)
#define SMF_ASSERT(x) ((void)0)
///FAIL_USING_EXCEPTIONS WILL CAUSE MathException to be raised on fail.
#elif defined(FAIL_USING_EXCEPTIONS)

#include <stdexcept>

#define SMF_ASSERT(x) \
		if (!(x)){ \
		string error(#x);  \
		SMF::Debug::debug(SMF::Debug::error,__FUNCTION__) <<"Assumption "  << #x << " failed! in file "<<  __FILE__ << " line  " << __LINE__ << endl; \
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);\
	}\


#elif defined(_MSC_VER)  //SMF_ASSERT will not generate exceptions, will only generate debug messages

#define SMF_ASSERT(x) (void)((!!(x)) || ( SMF::Debug::debug(SMF::Debug::error,__FUNCTION__) <<"Assumption "  << #x << " failed! in file "<<  __FILE__ << " line  " << __LINE__ << endl;

#elif defined(ANDROID)

#include <android/log.h>
#define SMF_ASSERT(x) do { if (!(x)) { __android_log_print(ANDROID_LOG_ERROR, "native-activity", "Assumption \"%s\" failed! in file %s, line %d!\n", #x, __FILE__, __LINE__); } } while(0)
#ifdef SMF_ASSERT
#undef SMF_ASSERT
#endif
#define SMF_ASSERT(x) do { if (!(x)) { __android_log_print(ANDROID_LOG_ERROR, "native-activity", "Assertion \"%s\" failed! in file %s, line %d!\n", #x, __FILE__, __LINE__); } } while(0)

#else // All other platforms

#define SMF_ASSERT(x) do { if (!(x)) { SMF:Debug::debug(SMF::Debug::error,__FUNCTION__) <<"Assumption "  << #x << " failed! in file "<<  __FILE__ << " line  " << __LINE__ << endl; } } while(0)

#endif

// If TEST_FOR_CORRECTNESS is defined, the function SMF_ASSERT() is enabled to test
// that all forms of optimizations inside the math library produce proper results.
#define TEST_FOR_CORRECTNESS


// Kill SMF_ASSERT() macros in OPTIMIZED_RELEASE builds.
#ifdef OPTIMIZED_RELEASE
#ifdef SMF_ASSERT
#undef SMF_ASSERT
#endif
#define SMF_ASSERT(x) ((void)0)
#define mathassert(x) ((void)0)
#endif







    #define assertmem( x, y )				SMF_ASSERT( _CrtIsValidPointer( x, y, true ) )
    #ifndef ASSERT
    #if defined(DEBUG) || defined(_DEBUG)
    #define ASSERT(test) if( !(test) ) { \
	fprintf(stderr,"\nASSERT(%s) FAILS, %s line %d\n",#test,__FILE__, __LINE__); exit(0);}
    #else
    #define ASSERT(test)
    #endif
    #endif
/**
To check this, cast the pointer to an integer of suitable size, take the modulus ALIGN, and check whether the result is zero.
 In code:
 bool is_aligned(void *p, int ALIGN)
{
    return (int)p % ALIGN == 0;
}
 To check for alignment of a power of 2 you can use:

((unsigned long)p & (ALIGN - 1)) == 0

This is simply a faster version of (p % ALIGN) == 0.

(If ALIGN is a constant your compiler will probably automatically use the faster version above.)

**/
#define assert_X_byte_aligned( ptr, x )		SMF_ASSERT( ( ((UINT_PTR)(ptr)) &  (x-1) ) == 0 )
#define assert_2_byte_aligned( ptr )		SMF_ASSERT( ( ((UINT_PTR)(ptr)) &  1 ) == 0 )
#define assert_4_byte_aligned( ptr )		SMF_ASSERT( ( ((UINT_PTR)(ptr)) &  3 ) == 0 )
#define assert_8_byte_aligned( ptr )		SMF_ASSERT( ( ((UINT_PTR)(ptr)) &  7 ) == 0 )
#define assert_16_byte_aligned( ptr )		SMF_ASSERT( ( ((UINT_PTR)(ptr)) & 15 ) == 0 )
#define assert_32_byte_aligned( ptr )		SMF_ASSERT( ( ((UINT_PTR)(ptr)) & 31 ) == 0 )
#define assert_64_byte_aligned( ptr )		SMF_ASSERT( ( ((UINT_PTR)(ptr)) & 63 ) == 0 )
#define assert_128_byte_aligned( ptr )		SMF_ASSERT( ( ((UINT_PTR)(ptr)) & 127 ) == 0 )
#define assert_aligned_to_type_size( ptr )	SMF_ASSERT( ( ((UINT_PTR)(ptr)) & ( sizeof( (ptr)[0] ) - 1 ) ) == 0 )


template<bool> struct compile_time_assert_failed;
template<> struct compile_time_assert_failed<true> {};
template<int x> struct compile_time_assert_test {};
#define compile_time_assert_join2( a, b )	a##b
#define compile_time_assert_join( a, b )	compile_time_assert_join2(a,b)
#define compile_time_assert( x )			typedef compile_time_assert_test<sizeof(compile_time_assert_failed<(bool)(x)>)> compile_time_assert_join(compile_time_assert_typedef_, __LINE__)

#define assert_sizeof( type, size )						compile_time_assert( sizeof( type ) == size )
#define assert_sizeof_8_byte_multiple( type )			compile_time_assert( ( sizeof( type ) &  7 ) == 0 )
#define assert_sizeof_16_byte_multiple( type )			compile_time_assert( ( sizeof( type ) & 15 ) == 0 )
#define assert_offsetof( type, field, offset )			compile_time_assert( offsetof( type, field ) == offset )
#define assert_offsetof_8_byte_multiple( type, field )	compile_time_assert( ( offsetof( type, field ) & 7 ) == 0 )
#define assert_offsetof_16_byte_multiple( type, field )	compile_time_assert( ( offsetof( type, field ) & 15 ) == 0 )

/////MEMORY HEAP DEFINITIONS///////////////
// uncoment next line to debug SMF heap memory utilization
//#define SMF_DEBUG_MEM

// uncomnt next line to redirect new and delete to SMF new and SMF delete
//#define SMF_REDIRECT_NEWDELETE

#endif //end WINDOWS













//////////////////////////////////////////////////////////
//  INLINE DEFINITIONS
/////////////////////////////////////////////////////////
#ifndef SMF_INLINE
#if defined(_MSC_VER)
#define SMF_INLINE __inline
#else //mingw
#define SMF_INLINE __inline__
#endif

#endif
#ifndef SMF_INLINE_FORCED

    #if defined(_MSC_VER)
        #define SMF_INLINE_FORCED __forceinline
    #else  //Mingw
        //#define SMF_INLINE_FORCED __attribute__((always_inline))
        #define SMF_INLINE_FORCED __inline__

    #endif
#endif

#ifndef SMF_INLINE_EXTERN
// lint complains that extern used with definition is a hazard, but it
// has the benefit (?) of making it illegal to take the address of the function
#ifdef _lint
#define SMF_INLINE_EXTERN			__forceinline
#else
#if defined (_MSC_VER)
#define SMF_INLINE_EXTERN			extern __forceinline
#else // MingW
#define SMF_INLINE_EXTERN			inline __attribute__((always_inline))
#endif
#endif



#endif

////////////////////////////////////////////////////////////
// Define portable import / export macros
////////////////////////////////////////////////////////////
#if defined(SMF_ON_WINDOWS)

    #ifdef SMF_DYNAMIC

        // Windows platforms
        #ifdef SGE_EXPORTS

            // From DLL side, we must export
            #define SMF_API __declspec(dllexport)

        #else

            // From client application side, we must import
            #define SMF_API __declspec(dllimport)

        #endif

        // For Visual C++ compilers, we also need to turn off this annoying C4251 warning.
        // You can read lots ot different things about it, but the point is the code will
        // just work fine, and so the simplest way to get rid of this warning is to disable it
        #ifdef _MSC_VER

            #pragma warning(disable : 4251)
/* Disable ridiculous warnings so that the code */
/* compiles cleanly using warning level 4.      */

/* nonstandard extension 'single line comment' was used */
#pragma warning(disable: 4001)

// indirection to slightly different base types
#pragma warning(disable: 4057)

// unreferenced formal parameter
#pragma warning(disable: 4100)

// named type definition in parentheses
#pragma warning(disable: 4115)

// nonstandard extension used : nameless struct/union
#pragma warning(disable: 4201)

// nonstandard extension used : benign typedef redefinition
#pragma warning(disable: 4209)

// nonstandard extension used : bit field types other than int
#pragma warning(disable: 4214)

// unreferenced inline function has been removed
#pragma warning(disable: 4514)

// C++ language change: to explicitly specialize class template '%s' use the following syntax:
#pragma warning(disable: 4663)

// Note: Creating precompiled header
#pragma warning(disable: 4699)

// Identifier was truncated to '255' characters in the debug information
#pragma warning(disable: 4786)

//--------------------------------------------------------------------------------------------------
// Do not tolerate the following L4 warnings:

// 'function' : different types for formal and actual parameter 'number'
#pragma warning( error: 4024 )

// formal parameter 'number' different from declaration
#pragma warning( error: 4028 )

// formal parameter 'number' has different type when promoted
#pragma warning( error: 4032 )

// 'identifier1' : 'operator' : different levels of indirection from 'identifier2'
#pragma warning( error: 4047 )

//enumerate 'identifier' in switch of enum 'identifier' is not explicitly handled by a case label
#pragma warning( error: 4061 )

// unreferenced formal parameter
//#pragma warning(disable: 4100)

// unreferenced local variable
//#pragma warning(disable: 4101)

// alignment of a member was sensitive to packing
#pragma warning( error: 4121 )

// logical operation on address of string constant
#pragma warning( error: 4130 )

// const object should be initialized
#pragma warning( error: 4132 )

// An attempt was made to subtract two pointers of different types. (?!)
#pragma warning( error: 4133 )

// local variable is initialized but not referenced
//#pragma warning( error: 4189 )

// nonstandard extension used : float long
#pragma warning( error: 4216 )

// nonstandard extension used : 'identifier' : cannot be initialized using address of automatic variable
#pragma warning( error: 4221 )

// cast truncates constant value
#pragma warning( error: 4310 )

// '==' : operator has no effect; did you intend '='?
#pragma warning( error: 4553 )

// local variable 'identifier' used without having been initialized
#pragma warning( error: 4700 )

// local variable 'identifier' may be used without having been initialized
//#pragma warning( error: 4701 )

// assignment within conditional expression
#pragma warning( error: 4706 )

// not all control paths return a value
#pragma warning( error: 4715 )



        #endif

    #else

        // No specific directive needed for static build
        #define SMF_API

    #endif

#else

    // Other platforms don't need to define anything
    #define SMF_API

#endif
////////////////////////////////////////////////////////////////////
//   outras definições vlidas para qualquer plataforma
///////////////////////////////////////////////////////////////////
#define	MAX_STRING_CHARS		1024		// max length of a string

#ifndef MARK_UNUSED
/// If a variable is labelled with this directive, the compiler should not emit a warning even if it is unused in the code.
#define MARK_UNUSED(x) ((void)x)
#endif

#ifdef __GNUC__
/// If a variable or a function definition is labelled with this directive, the compiler should not emit a warning even if it is unused
/// in the code.
#define DONT_WARN_UNUSED __attribute__((unused))
#else
#define DONT_WARN_UNUSED
#endif
//to teste macros
#define SMF_CONTAINS(status, x)  (status & x) != 0

#endif
