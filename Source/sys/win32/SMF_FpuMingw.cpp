/*
  SMF - Salvathor Math Fabric  (https://sourceforge.net/projects/smfabric/)
  Copyright (C) 2010-2011 Salvatore Giannotta Filho <a_materasu@hotmail.com>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

  */


#include "sys/SMF_System.h"
#if defined (WIN32)
#include <string>
//#include  <io.h>
#include  <stdio.h>
#include  <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdint.h>
#include "sys/SMF_System.h"
#include "util/SMF_Debug.h"
#include "util/SMF_StringUtils.h"
//#define _GNU_SOURCE         /* See feature_test_macros(7) */
#include <fenv.h>
#undef __STRICT_ANSI__ // _controlfp is a non-standard function documented in MSDN
#include <float.h>

#include <math.h>
#include <xmmintrin.h>

namespace SMF{

namespace System{
using namespace std;


/*
===============================================================================

	FPU - Float Point Unit
//http://www.website.masmforum.com/tutorials/fptute/fpuchap3.htm
===============================================================================
*/
//Exceptions
// control word
/*
 * Control Word

The Control Word 16-bit register is used by the programmer to select between the various modes of computation available from the FPU,
and to define which exceptions should be handled by the FPU or by an exception handler written by the programmer.
The Control Word is divided into several bit fields as depicted in the following Fig.1.2.
Bits 6,7,13,14,15 are reserved and not used.

Bits
EXCEPTION-FLAG MASKS 0-5
0 -> IM  Invalid Operation
1 -> DM  Denormalized Operand
2 -> ZM  Zero Divide
3 -> OM  Overflow
4 -> UM  Underflow
5 -> PM  Precision
PRECISION CONTROL FIELD
8,9 PC ->  Single Precision   (24-bit) = $00B
           Double Precision   (53-bit) = $10B
           Extended Precision (64-bit) = $11B
           Reserved                    = $01B

ROUNDING MODE
10,11 RC ->  Round to nearest even      = $00B
             Round down toward infinity = $01B
             Round up toward infinity   = $10B
             Round toward zero (trunc)  = $11B

INFINITY CONTROL
12  X ->     Infinity control. Used for compatability with 287 FPU.
             Projective = 0
             Affine     = 1


*/

#ifdef __ASSEMBLER__
#define _Const_(x)      $##x
#else
#define _Const_(x)      x
#endif

#define FPU_CW_RC           _Const_(0x0C00) /* rounding control */
#define FPU_CW_PC           _Const_(0x0300) /* precision control */
#define FPU_CW_PM           _Const_(0x0020) /* precision mask */
#define FPU_CW_UM           _Const_(0x0010) /* underflow mask */
#define FPU_CW_OM           _Const_(0x0008) /* overflow mask */
#define FPU_CW_ZM           _Const_(0x0004) /* divide by zero mask */
#define FPU_CW_DM           _Const_(0x0002) /* denormalized operand mask */
#define FPU_CW_IM           _Const_(0x0001) /* invalid operation mask */
#define FPU_CW_EXM          _Const_(0x007f) /* all masks */


//  =   STATUS WORD =======
/*Status word

The Status Word 16-bit register indicates the general condition of the FPU.
Its content may change after each instruction is completed.
Part of it cannot be changed directly by the programmer.
It can, however, be accessed indirectly at any time to inspect its content.
*/

#define FPU_SW_B            Const__(0x8000) /* backward compatibility (=ES) */
#define FPU_SW_C3           Const__(0x4000) /* condition bit 3 */
#define FPU_SW_TOP          Const__(0x3800) /* top of stack */
#define FPU_SW_TOPS         Const__(11)     /* shift for top of stack bits */
#define FPU_SW_C2           Const__(0x0400) /* condition bit 2 */
#define FPU_SW_C1           Const__(0x0200) /* condition bit 1 */
#define FPU_SW_C0           Const__(0x0100) /* condition bit 0 */
#define FPU_SW_ES           Const__(0x0080) /* exception summary */
#define FPU_SW_SF           Const__(0x0040) /* stack fault */
#define FPU_SW_PE           Const__(0x0020) /* loss of precision */
#define FPU_SW_UE           Const__(0x0010) /* underflow */
#define FPU_SW_OE           Const__(0x0008) /* overflow */
#define FPU_SW_ZE           Const__(0x0004) /* divide by zero */
#define FPU_SW_DE           Const__(0x0002) /* denormalized operand */
#define FPU_SW_IE           Const__(0x0001) /* invalid operation */

/*
 * Tag Word

The Tag Word 16-bit register is managed by the FPU to maintain some information on the content of each of its 80-bit registers.

The Tag Word is divided into 8 fields of 2 bits each as depicted in the following Fig.1.4.
*/
#define FPU_TW_TAG7	   Const__(0xC000)
#define FPU_TW_TAG6	   Const__(0x3000)
#define FPU_TW_TAG5	   Const__(0x0C00)
#define FPU_TW_TAG4	   Const__(0x0300)
#define FPU_TW_TAG3	   Const__(0x00C0)
#define FPU_TW_TAG2	   Const__(0x0030)
#define FPU_TW_TAG1	   Const__(0x000C)
#define FPU_TW_TAG0	   Const__(0x0003)

// Exceptions flags
#define FPU_BUSY        Const_(0x8000)   /* FPU busy bit (8087 compatibility) */
#define FPU_EX_ERRORSUMARY Const_(0x0080)   /* Error summary status */
/* Special exceptions: */
#define EFPU_X_INTERNAL     Const_(0x8000)  /* Internal error in wm-FPU-emu */
#define FPU_EX_STACKOVER    Const_(0x0041|SW_C1)    /* stack overflow */
#define FPU_EX_STACKUNDER   Const_(0x0041)  /* stack underflow */
// FLAGS
#define FPU_EX_INVALID		Const_(0x0001)  /* invalid operation */
#define FPU_EX_DENORM		Const_(0x0002)  /* denormalized operand */
#define FPU_EX_DIVZERO		Const_(0x0004)  /* divide by zero */
#define FPU_EX_OVERFLOW		Const_(0x0008)  /* overflow */
#define FPU_EX_UNDERFLOW	Const_(0x0010)  /* underflow */
#define FPU_EX_PRECISION	Const_(0x0020)  /* loss of precision */ // INEXACT

#define PARAM1(x)  8(x)

/* rOUNDING cONTROL
 */

#define FPU_ROU_TONEAREST FE_TONEAREST
#define FPU_ROU_UPWARD FE_UPWARD
#define FPU_ROU_DOWNWARD FE_DOWNWARD
#define FPU_ROU_TOWARDZERO FE_TOWARDZERO
/*
Precision Control
* 11 - round to extended precision
* 10 - round to double precision
* 00 - round to single precision
*/
/* precision control */
#define FPU_EXTENDED 0x300	//Single precision
#define FPU_DOUBLE   0x200     // Double precision
#define FPU_SINGLE   0x0       // Double-extended precision



typedef struct bitFlag_s {
	char *	name;
	int	bit;
} bitFlag_t;
#if defined(MASM_INTEL)
static sf_u8 fpuState[128];
static void *statePtr=fpuState;

#else
static fenv_t fpuState;
static fenv_t *statePtr=NULL;
#endif
//static statePtr = fpuState;

static char fpuString[2048];

static bitFlag_t controlWordFlags[] = {
	{ "Invalid operation", 0x01 },
	{ "Denormalized operand", 0x02 },
	{ "Divide-by-zero", 0x04 },
	{ "Numeric overflow", 0x08 },
	{ "Numeric underflow", 0x10 },
	{ "Inexact result (precision)", 0x20 },
	{ "Infinity control", 12 },
	{ "", 0 }
};
static char *precisionControlField[] = {
	"Single Precision (24-bits)",
	"Reserved",
	"Double Precision (53-bits)",
	"Double Extended Precision (64-bits)"
};

static char *roundingControlField[] = {
	"Round to nearest",
	"Round down",
	"Round up",
	"Round toward zero"
};
static bitFlag_t statusWordFlags[] = {
	{ "Invalid operation", 0 },
	{ "Denormalized operand", 1 },
	{ "Divide-by-zero", 2 },
	{ "Numeric overflow", 3 },
	{ "Numeric underflow", 4 },
	{ "Inexact result (precision)", 5 },
	{ "Stack fault", 6 },
	{ "Error summary status", 7 },
	{ "FPU busy", 15 },
	{ "", 0 }
};


static sf_u32 invert(sf_u32 x)
{
   sf_u32 y;
   y = x ^ 0xffffffff ;    // Do XOR to invert the bits , this will work for 32 bit compiler because 0xffffffff represents a 32-bit number with all bits set to 1.
   return (y);
   }
/*
===============
Sys_FPU_PrintStateFlags
===============
*/
int CFPU::fpu_PrintStateFlags( char *ptr, int ctrl, int stat, int tags, int inof, int inse, int opof, int opse ) {
	int i, length = 0;

	length += sprintf( ptr+length,	"CTRL = %08x\n"
					"STAT = %08x\n"
					"TAGS = %08x\n"
					"INOF = %08x\n"
					"INSE = %08x\n"
					"OPOF = %08x\n"
					"OPSE = %08x\n"
					"\n",
					ctrl, stat, tags, inof, inse, opof, opse );

	length += sprintf( ptr+length, "Control Word:\n" );
	for ( i = 0; controlWordFlags[i].name[0]; i++ ) {
		length += sprintf( ptr+length, "  %-30s = %s\n", controlWordFlags[i].name, ( ctrl & ( 1 << controlWordFlags[i].bit ) ) ? "true" : "false" );
	}
	length += sprintf( ptr+length, "  %-30s = %s\n", "Precision control", precisionControlField[(ctrl>>8)&3] );
	length += sprintf( ptr+length, "  %-30s = %s\n", "Rounding control", roundingControlField[(ctrl>>10)&3] );

	length += sprintf( ptr+length, "Status Word:\n" );
	for ( i = 0; statusWordFlags[i].name[0]; i++ ) {
		ptr += sprintf( ptr+length, "  %-30s = %s\n", statusWordFlags[i].name, ( stat & ( 1 << statusWordFlags[i].bit ) ) ? "true" : "false" );
	}
	length += sprintf( ptr+length, "  %-30s = %d%d%d%d\n", "Condition code", (stat>>8)&1, (stat>>9)&1, (stat>>10)&1, (stat>>14)&1 );
	length += sprintf( ptr+length, "  %-30s = %d\n", "Top of stack pointer", (stat>>11)&7 );

	return length;
}

/*
===============
Sys_FPU_StackIsEmpty
===============
*/
bool CFPU::fpu_StackIsEmpty( void ) {




	sf_u32 y;
	#if defined(MASM_INTEL)
__asm__ (
		"mov			eax, %0\n"
		"fnstenv		[eax]\n"
		"mov			eax, [eax+8]\n"
		"xor			eax, 0xFFFFFFFF\n"
		"and			eax, 0x0000FFFF\n"
		"jz			1\n"
	::"r"(statePtr)
    :);

    #else
	statePtr = &fpuState;
	fegetenv(statePtr);
	//status word = 16 bits
	asm volatile("movl (%0), %%eax\n\t"
		      : "+S" (statePtr));
	asm volatile("movl 8(%eax), %eax\n\t" //status word
		     "xor $0xFFFFFFFF , %eax\n\t"
		     "and $0x0000FFFF , %eax\n\t"
		     "jz 1f\n");
	#endif
	return false;
	asm volatile("1:");
	return true;

}

/*
===============
Sys_FPU_ClearStack
===============
*/
void CFPU::fpu_ClearStack( void ) {
#if defined(MASM_INTEL)
	__asm__ (
		"mov			eax, %0\n"
		"fnstenv		[eax]\n"
		"mov			eax, [eax+8]\n"
		"xor			eax, 0xFFFFFFFF\n"
		"mov			edx, (3<<14)\n"
	"emptyStack%=:\n"
		"mov			ecx, eax\n"
		"and			ecx, edx\n"
		"jz			done1%=\n"
		"fstp		st\n"
		"shr			edx, 2\n"
		"jmp			emptyStack%=\n"
	"done1%=:\n"
	::"r"(statePtr)
    :);
#else

	statePtr = &fpuState;
	fegetenv(statePtr);
/*  TODO
		asm volatile("movl (%0), %%eax\n\t"
		      : "+S" (statePtr));

		asm volatile("movl 8(%eax), %eax\n\t" //status word
		      "xor $0xFFFFFFFF , %eax\n\t"
		      "movl $0xC000 , %edx\n\t"
		      "1: movl %eax , %ecx\n\t"
		      "and %edx , % ecx \n\t"
		      "jz 2f\n\t"
		      "fst  st\n\t"    //although this instruction is allowed, no data register is changed
					//an exception may still be detected due to the value in ST(0)
		      "shr $2 , %edx\n\t"
		      "jz 1f\n\t"

		      "2:\n\t");
		      */
#if 0
	return true;

	__asm__ __volatile__(
		mov		eax, statePtr
		fnstenv		[eax]
		mov		eax, [eax+8] //status word
		xor		eax, 0xFFFFFFFF
		mov		edx, (3<<14)
	emptyStack:
		mov		ecx, eax
		and		ecx, edx
		jz		done
		fstp		st
		shr		edx, 2
		jmp		emptyStack
	done:
		)
#endif
#endif
}

//get the content of FPU registers without changing them
// return the number of values on the stack
static int getFPURegisters(double *var){

  sf_u16 status;
  int numValues=0;
unsigned char *Src1;
 unsigned char *Dest;
  unsigned int SrcLength;
   unsigned char C;
   	unsigned char *Mask;

#if defined(MASM_INTEL)

#else
	asm volatile
		("pusha		     \n\t"
		/* ** Duplicate C in 8 bytes of MM1 ** */
		"mov           %3, %%al \n\t"	/* load C into AL */
		"mov         %%al, %%ah \n\t"	/* copy AL into AH */
		"mov         %%ax, %%bx \n\t"	/* copy AX into BX */
		"shl         $16, %%eax \n\t"	/* shift 2 bytes of EAX left */
		"mov         %%bx, %%ax \n\t"	/* copy BX into AX */
		"movd      %%eax, %%mm1 \n\t"	/* copy EAX into MM1 */
		"movd      %%eax, %%mm2 \n\t"	/* copy EAX into MM2 */
		"punpckldq %%mm2, %%mm1 \n\t"	/* fill higher bytes of MM1 with C */
		"movl         %4, %%edx \n\t"	/* load Mask address into edx */
		"movq    (%%edx), %%mm0 \n\t"	/* load Mask into mm0 */
		"mov          %1, %%eax \n\t"	/* load Src1 address into eax */
		"mov          %0, %%edi \n\t"	/* load Dest address into edi */
		"mov          %2, %%ecx \n\t"	/* load loop counter (SIZE) into ecx */
		"shr          $3, %%ecx \n\t"	/* counter/8 (MMX loads 8 bytes at a time) */
		".align 16              \n\t"	/* 16 byte alignment of the loop entry */
		"1:                     \n\t"
		"movq    (%%eax), %%mm2 \n\t"	/* load 8 bytes from Src1 into MM2 */
		"psrlw        $1, %%mm2 \n\t"	/* shift 4 WORDS of MM2 1 bit to the right */
		/*    "pand      %%mm0, %%mm2 \n\t"    // apply Mask to 8 BYTES of MM2 */
		".byte     0x0f, 0xdb, 0xd0 \n\t"
		"paddusb   %%mm1, %%mm2 \n\t"	/* MM2=SrcDest+C (add 8 bytes with saturation) */
		"movq    %%mm2, (%%edi) \n\t"	/* store result in Dest */
		"add          $8, %%eax \n\t"	/* increase Src1 register pointer by 8 */
		"add          $8, %%edi \n\t"	/* increase Dest register pointer by 8 */
		"dec              %%ecx \n\t"	/* decrease loop counter */
		"jnz                  1b \n\t"	/* check loop termination, proceed if required */
		"emms                   \n\t"	/* exit MMX state */
		"popa                   \n\t":"=m" (Dest)	/* %0 */
		:"m"(Src1),		/* %1 */
		"m"(SrcLength),		/* %2 */
		"m"(C),			/* %3 */
		"m"(Mask)			/* %4 */
		);

/*
  asm volatile ("fnstsw %0" : "=a" (status));
  asm volatile ("movl var, %edi\n\t");
  asm volatile ("mov status, %esi\n\t");
  asm volatile ("xor $0xFFFFFFFF, %esi\n\t"
		  "movl $0xC000 , %edx\n\t"
		  "mov %esi, %ecx\n\t"
		  "and %edx,%ecx\n\t"
		  "jz done\n\t"

  );
  asm volatile ("fstpl %0;"     			//The fst and fstp instructions copy the value on the top of the floating point register stack to another floating point register or to a 32, 64, or 80 bit memory variable
		"inc %%eax\n\t"   		//The fst and fstp instructions will set the stack exception bit if a stack underflow occurs (attempting to store a value from an empty register stack)
		"shr $2, %%edx\n\t"
		"movl %%esi, %%ecx\n\t"
		"and  %%edx,%%ecx\n\t" 		//The B field (bit 15) indicates if the FPU is busy (B=1) while executing an instruction, or is idle (B=0).
						//The C3 (bit 14)  fields contain the condition codes following the execution of some instructions such as comparisons
						// zero = if bit 14 and 15 are 1 originally (before invertion). ie . The FPU is Busy and C3 is
		"jz done\n\t": "=m"(var[0]):);

  asm volatile(	"fxch %%st(1)\n\t" 		//The fxch instruction exchanges the value on the top of stack with st(1)
		"fstpl %0;"
		"inc %%eax\n\t"
		"fxch %%st(1)\n\t"
		"shr $2, %%edx\n\t"		// mudaa a mascara para os bits
		"movl %%esi, %%ecx\n\t"
		"and  %%edx,%%ecx\n\t"
		"jz done\n\t": "=m"(var[1]):);

  asm volatile(	"fxch %%st(2)\n\t" 		//The fxch instruction exchanges the value on the top of stack with st(2)
		"fstpl %0;"
		"inc %%eax\n\t"
		"fxch %%st(2)\n\t"
		"shr $2, %%edx\n\t"		// mudaa a mascara para os bits
		"movl %%esi, %%ecx\n\t"
		"and  %%edx,%%ecx\n\t"
		"jz done\n\t": "=m"(var[2]):);
  asm volatile(	"fxch %%st(3)\n\t" 		//The fxch instruction exchanges the value on the top of stack with st(2)
		"fstpl %0;"
		"inc %%eax\n\t"
		"fxch %%st(3)\n\t"
		"shr $2, %%edx\n\t"		// mudaa a mascara para os bits
		"movl %%esi, %%ecx\n\t"
		"and  %%edx,%%ecx\n\t"
		"jz done\n\t": "=m"(var[3]):);
  asm volatile(	"fxch %%st(4)\n\t" 		//The fxch instruction exchanges the value on the top of stack with st(2)
		"fstpl %0;"
		"inc %%eax\n\t"
		"fxch %%st(4)\n\t"
		"shr $2, %%edx\n\t"		// mudaa a mascara para os bits
		"movl %%esi, %%ecx\n\t"
		"and  %%edx,%%ecx\n\t"
		"jz done\n\t": "=m"(var[4]):);
  asm volatile(	"fxch %%st(5)\n\t" 		//The fxch instruction exchanges the value on the top of stack with st(2)
		"fstpl %0;"
		"inc %%eax\n\t"
		"fxch %%st(5)\n\t"
		"shr $2, %%edx\n\t"		// mudaa a mascara para os bits
		"movl %%esi, %%ecx\n\t"
		"and  %%edx,%%ecx\n\t"
		"jz done\n\t": "=m"(var[5]):);
  asm volatile(	"fxch %%st(6)\n\t" 		//The fxch instruction exchanges the value on the top of stack with st(2)
		"fstpl %0;"
		"inc %%eax\n\t"
		"fxch %%st(6)\n\t"
		"shr $2, %%edx\n\t"		// mudaa a mascara para os bits
		"movl %%esi, %%ecx\n\t"
		"and  %%edx,%%ecx\n\t"
		"jz done\n\t": "=m"(var[6]):);
    asm volatile("fxch %%st(7)\n\t" 		//The fxch instruction exchanges the value on the top of stack with st(2)
		"fstpl %0;"
		"inc %%eax\n\t"
		"fxch %%st(7)\n\t": "=m"(var[7]):);
    asm volatile ("done:\n\t"
		"movl	%eax, numValues\n\t");
*/
#endif
}

/*
===============
fpu_GetState

  gets the FPU state without changing the state
===============
*/
const char *CFPU::fpu_GetState( void ) {
	double fpuStack[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double *fpuStackPtr = fpuStack;
	int i, numValues;
	char *ptr;

#if defined(MASM_INTEL)

	__asm__ __volatile__(
		"mov			esi, %0\n"
		"mov			edi, %1\n"
		"fnstenv		[esi]\n"
		"mov			esi, [esi+8]\n" 		//status word
		"xor			esi, 0xFFFFFFFF	\n"	// inverte status word
		"mov			edx, (3<<14)\n"           // mascara
		"xor			eax, eax\n"        	// zera eax ; contador de nro de valores
		"mov			ecx, esi\n"
		"and			ecx, edx\n"
		"jz			done\n"
		"fst			qword ptr [edi+0]\n"   	//The fst and fstp instructions copy the value on the top of the floating point register stack to another floating point register or to a 32, 64, or 80 bit memory variable
		"inc			eax\n"			//The fst and fstp instructions will set the stack exception bit if a stack underflow occurs (attempting to store a value from an empty register stack)
		"shr			edx, 2\n"
		"mov			ecx, esi\n"
		"and			ecx, edx \n"           //The B field (bit 15) indicates if the FPU is busy (B=1) while executing an instruction, or is idle (B=0).
							    //The C3 (bit 14)  fields contain the condition codes following the execution of some instructions such as comparisons
							    // zero = if bit 14 and 15 are 1 originally (before invertion). ie . The FPU is Busy and C3 is
		"jz			done\n"
		"fxch		st(1)\n"    //The fxch instruction exchanges the value on the top of stack with st(1)
		"fst			qword ptr [edi+8]\n"
		"inc			eax\n"
		"fxch		st(1)\n"
		"shr			edx, 2\n"     	    // mudaa a mascara para os bits
		"mov			ecx, esi\n"
		"and			ecx, edx\n"
		"jz			done\n"
		"fxch		st(2)\n"
		"fst			qword ptr [edi+16]\n"
		"inc			eax\n"
		"fxch		st(2)\n"
		"shr			edx, 2\n"
		"mov			ecx, esi\n"
		"and			ecx, edx\n"
		"jz			done\n"
		"fxch		st(3)\n"
		"fst			qword ptr [edi+24]\n"
		"inc			eax\n"
		"fxch		st(3)\n"
		"shr			edx, 2\n"
		"mov			ecx, esi\n"
		"and			ecx, edx\n"
		"jz			done\n"
		"fxch		st(4)\n"
		"fst			qword ptr [edi+32]\n"
		"inc			eax\n"
		"fxch		st(4)\n"
		"shr			edx, 2\n"
		"mov			ecx, esi\n"
		"and			ecx, edx\n"
		"jz			done\n"
		"fxch		st(5)\n"
		"fst			qword ptr [edi+40]\n"
		"inc			eax\n"
		"fxch		st(5)\n"
		"shr			edx, 2\n"
		"mov			ecx, esi\n"
		"and			ecx, edx\n"
		"jz			done\n"
		"fxch		st(6)\n"
		"fst			qword ptr [edi+48]\n"
		"inc			eax\n"
		"fxch		st(6)\n"
		"shr			edx, 2\n"
		"mov			ecx, esi\n"
		"and			ecx, edx\n"
		"jz			done\n"
		"fxch		st(7)\n"
		"fst			qword ptr [edi+56]\n"
		"inc			eax\n"
		"fxch		st(7)\n"
	"done:\n"
		"mov		%2	, eax\n"
	::"r"(statePtr),"r"(fpuStackPtr),"m"(numValues)
    :);

	int ctrl = *(int *)&fpuState[0];
	int stat = *(int *)&fpuState[4];
	int tags = *(int *)&fpuState[8];
	int inof = *(int *)&fpuState[12];
	int inse = *(int *)&fpuState[16];
	int opof = *(int *)&fpuState[20];   //operand address
	int opse = *(int *)&fpuState[24];   //data selector


#else
	statePtr = &fpuState;

	fegetenv(statePtr);

	numValues = getFPURegisters(fpuStack);

#endif
	ptr = fpuString;
	ptr += sprintf( ptr,"FPU State:\n"
						"num values on stack = %d\n", numValues );
	for ( i = 0; i < 8; i++ ) {
		ptr += sprintf( ptr, "ST%d = %1.10e\n#define __SSE__", i, fpuStack[i] );
	}
#if defined(MASM_INTEL)

#else
	Sys_FPU_PrintStateFlags( ptr, fpuState.__control_word, fpuState.__status_word, fpuState.__tag_word, fpuState.__ip_offset, fpuState.__ip_selector, fpuState.__opcode, fpuState.__data_selector );
#endif
	return fpuString;
}

/*
===============
fpu_EnableExceptions

===============
*/



void CFPU::fpu_EnableExceptions( fpuExceptions_t exceptions ) {
//_controlfp

_clearfp();
   unsigned unused_current_word = 0;
   _controlfp_s(&unused_current_word, 0, exceptions);  // _controlfp_s is the secure version of _controlfp

//  if(feenableexcept(exceptions) == -1)
 //   printf("erro");
}
/*
===============
fpu_DisableExceptions
===============
*/
void CFPU::fpu_DisableExceptions( int excepts ) {

     // This code is what is suggested in Solaris docs
// I'm guessing that was a cruel hoax

 int e = 0;
 if((excepts & FE_DIVBYZERO) != 0)
 {
 e &= ~FE_DIVBYZERO;
 }
else if((excepts & FE_INVALID) != 0)
{
 e &= ~FE_INVALID;
}
else if((excepts & FE_OVERFLOW) != 0)
 {
 e &= ~FE_OVERFLOW;
}
 else if((excepts & FE_UNDERFLOW) != 0)
{
 e &= ~FE_UNDERFLOW;
 }
else if((excepts & FE_DENORMAL) != 0)
{
 e &= ~FE_DENORMAL;
 }
 else if((excepts & FE_INEXACT) != 0)
{
 e &= ~FE_INEXACT;
 }
  unsigned unused_current_word = 0;
   _controlfp_s(&unused_current_word, 0, e);  // _controlfp_s is the secure version of _controlfp

}
/*
===============
Sys_FPU_SetPrecision
===============
/* precision control
Param int precision
0 = Double extended Precision
1 = Undefined
2 = Double Precision
3 = Extended Precision
*/
void CFPU::fpu_SetPrecision( fpuPrecision_t precision ) {
	//short precisionBitTable[4] = { 0, 1, 3, 0 };
	//short precisionBits = precisionBitTable[precision & 3] << 8;
	short precisionBits = (precision & 3) << 8;
	short precisionMask = ~(FPU_CW_PC);//~( ( 1 << 9 ) | ( 1 << 8 ) ); //
	sf_u16 _cw;
	/*
	1 - get control word
	2 - zera os bits de precision
	3 - adiciona a nova precision na control word com .OR.
	4 - chama o comando fldcw para carregar a nova control word
	*/
	asm volatile ("fnstcw %0;": "=m" (_cw));
	_cw = (_cw && precisionMask);
	_cw = (_cw || precisionBits);
	asm volatile ("fldcw %0" : : "m" (*&_cw));
	/*
	__asm__ __volatile__(
		mov		eax, statePtr
		mov		cx, precisionBits
		fnstcw		word ptr [eax]  //Stores the current value of the FPU control word at the specified destination in memory
		mov		bx, word ptr [eax]
		and		bx, precisionMask //zera apenas os bits de precision na control word
		or		bx, cx
		mov		word ptr [eax], bx
		fldcw		word ptr [eax]
		)*/
}

void	CFPU::fpu_SetRounding( fpuRounding_t rounding ){
  if(fesetround(rounding)!=0)
  {
    printf("erro");
  }
}

int16_t	CFPU::fpu_getRounding( int16_t rounding ){
  return fegetround();
}
#if 0

/* Special i386 GCC implementation */
#if defined(__i386__) && defined(__GNUC__) && !defined(__BEOS__) && !defined(__FORCE_SSE2__)
#  define FPU_CONTROL
/* both GCC and MSVC are kinda stupid about rounding/casting to int.
   Because of encapsulation constraints (GCC can't see inside the asm
   block and so we end up doing stupid things like a store/load that
   is collectively a noop), we do it this way */

/* we must set up the fpu before this works!! */

typedef int16_t fpu_control;

static inline void fpu_setround(fpu_control *fpu){
  int16_t ret;
  int16_t temp;
  __asm__ __volatile__("fnstcw %0\n\t"
	  "movw %0,%%dx\n\t"
	  "orw $62463,%%dx\n\t"
	  "movw %%dx,%1\n\t"
	  "fldcw %1\n\t":"=m"(ret):"m"(temp): "dx");
  *fpu=ret;
}

static inline void fpu_restore(fpu_control fpu){
  __asm__ __volatile__("fldcw %0":: "m"(fpu));
}

/* assumes the FPU is in round mode! */
static inline int ftoi(double f){  /* yes, double!  Otherwise,
                                             we get extra fst/fld to
                                             truncate precision */
  int i;
  __asm__ __volatile__("fistl %0": "=m"(i) : "t"(f));
  return(i);
}
#endif /* Special i386 GCC implementation */


/* MSVC inline assembly. 32 bit only; inline ASM isn't implemented in the
 * 64 bit compiler */
#if defined(_MSC_VER) && !defined(_WIN64)
#  define FPU_CONTROL

typedef int16_t fpu_control;

static __inline int ftoi(double f){
	int i;
	__asm{
		fld f
		fistp i
	}
	return i;
}



/*
================
Sys_FPU_SetRounding
================
*/
__inline void Sys_FPU_SetRounding( int16_t rounding ) {

	int16_t roundingBitTable[4] = { 0, 1, 2, 3 };
	int16_t roundingBits = roundingBitTable[rounding & 3] << 10;
	int16_t roundingMask = ~( ( 1 << 11 ) | ( 1 << 10 ) );

	__asm__ __volatile__ (
		mov			eax, statePtr
		mov			cx, roundingBits
		fnstcw		word ptr [eax]
		mov			bx, word ptr [eax]
		and			bx, roundingMask
		or			bx, cx
		mov			word ptr [eax], bx
		fldcw		word ptr [eax]
		)
}


__inline void FPU_restoreRouding(fpu_control fpu){
}

#endif /* Special MSVC 32 bit implementation */

/* Optimized code path for x86_64 builds. Uses SSE2 intrinsics. This can be
   done safely because all x86_64 CPUs supports SSE2. */
#if (!defined(__FORCE_SSE2__)) || (defined(_MSC_VER) && defined(_WIN64)) || (defined(__GNUC__) && defined (__x86_64__))
#pragma warning "using sse2 for ftoi"
#define FPU_CONTROL

typedef int16_t fpu_control;

#include <emmintrin.h>
static __inline int ftoi(double f){
        return _mm_cvtsd_si32(_mm_load_sd(&f));
}

static __inline void fpu_setround(fpu_control *fpu){
}

static __inline void fpu_restore(fpu_control fpu){
}

#endif /* Special MSVC x64 implementation */


/* If no special implementation was found for the current compiler / platform,
   use the default implementation here: */
#ifndef FPU_CONTROL

typedef int fpu_control;

static int ftoi(double f){
        /* Note: MSVC and GCC (at least on some systems) round towards zero, thus,
           the floor() call is required to ensure correct roudning of
           negative numbers */
        return (int)floor(f+.5);
}

/* We don't have special code for this compiler/arch, so do it the slow way */
#  define fpu_setround(fpu_control) {}
#  define fpu_restore(fpu_control) {}

#endif /* default implementation */

#endif
#define X87FLAGBITS         6
#define DAZ_BIT            6
#define FTZ_BIT            15
#define DENORMAL_EXCEPTION_MASK   8
#define UNDERFLOW_EXCEPTION_MASK   11


// Esté código só compila comsuporte SSE. Caso contrário haverá erro

static void set_mxcsr_on(int bit_num){
  __m128    state[32];

  sf_u32   x;

  int error; // variable to get error form fxsave

CFPU::fpu_getFXSave ((void *) &state);
 memcpy( (void*)&x, (char*)state+24, 4);
x |= (1 << bit_num);
//asm volatile("ldmxcsr %0" : "=m" (x));
setmxcsr(x)  //LDMXCSR--Load Streaming SIMD Extension Control/Status

 }

static void set_mxcsr_off(int bit_num){


__m128    state[32];
sf_u32   x;
int error;  // variable to get error form fxsave
CFPU::fpu_getFXSave ((void *) &state);
memcpy( (void*)&x, (char*)state+24, 4);
x &= ~(1 << bit_num);
setmxcsr(x) //LDMXCSR--Load Streaming SIMD Extension Control/Status*/

 }

static void clear_flags(){
__m128    state[32];
CFPU::fpu_getFXSave ((void *) &state);
sf_u32  x;
int error;  // variable to get error form fxsave
//saveFX(x,error)
 memcpy( (void*)&x, (char*)state+24, 4);
x = x >> X87FLAGBITS;
x = x << X87FLAGBITS;
setmxcsr(x)

 }



// UNDERFLOWS
void CFPU::fpu_enableFTZ(){
  // UNDERFLOWS
set_mxcsr_on(FTZ_BIT);
set_mxcsr_off(UNDERFLOW_EXCEPTION_MASK);
}
void CFPU::fpu_disableFTZ(){
clear_flags();
}

// DENORMALS

void	CFPU::fpu_enableDAZ() {
set_mxcsr_off(DAZ_BIT);
set_mxcsr_on(DENORMAL_EXCEPTION_MASK);


}
void    CFPU::fpu_disableDAZ(){
clear_flags();
}

/*
================
Sys_FPU_SetDAZ
// http://software.intel.com/en-us/articles/x87-and-sse-floating-point-assists-in-ia-32-flush-to-zero-ftz-and-denormals-are-zero-daz
//Denormals-Are-Zero (DAZ)
DAZ is very similar to FTZ in many ways. DAZ mode is a method of bypassing IEEE 754 methods of dealing with denormal floating-point numbers. As mentioned, this mode is less precise, but much faster and is typically used in applications like streaming media when minute differences in quality are essentially undetectable. Two conditions must be met in order for DAZ processing to occur:

    The DAZ bit (bit 6) in the MXCSR register must be masked (value = 1).
    The processor must support DAZ mode. Initial steppings of Pentium® 4 processors did not support DAZ.
    Information about detection of DAZ can be found in the Intel® 64 and IA-32 Architectures Developer's Manual: Vol. 2A in the References section.


================
*/
void CFPU::fpu_SetDAZ( bool enable ) {
	sf_u32 dwData;

getmxcsr(dwData);
dwData = dwData && ~(1<<6);
dwData = dwData || ((enable && 1)<<6);
setmxcsr(dwData)

/*
	asm (
		"movzbl	 %1, %%ecx\n"
		"and		 1,%%ecx\n"
		"shl		 $6, %%ecx\n"   //DAZ bit =6
		"STMXCSR	%0\n"
		"mov		%0, %%eax\n"
		"and		 $~(1<<6),%%eax\n"	// clear DAX bit
		"or		    %%ecx, %%eax\n"		// set the DAZ bit
		"mov		%%eax, %0\n"
		"LDMXCSR	%0\n"
	::"m"(dwData),"m"(enable)
    :);

	__asm__ __volatile__ (
		movzx	ecx, byte ptr enable
		and		ecx, 1
		shl		ecx, 6         //DAZ bit =6
		STMXCSR	dword ptr dwData       //Stores the contents of the MXCSR control and status register to the destination operand. The destination operand is a 32-bit memory location. The reserved bits in the MXCSR register are stored as 0s.
		mov		eax, dwData
		and		eax, ~(1<<6)	// clear DAX bit
		or		eax, ecx	// set the DAZ bit
		mov		dwData, eax
		LDMXCSR	dword ptr dwData
		)*/
}

/*
================
Sys_FPU_SetFTZ
//http://software.intel.com/en-us/articles/x87-and-sse-floating-point-assists-in-ia-32-flush-to-zero-ftz-and-denormals-are-zero-daz
// Flush-To-Zero (FTZ)
================
*/
void fpu_SetFTZ( bool enable ) {

	sf_u32 dwData;
getmxcsr(dwData);
dwData = dwData && ~(1<<15);
dwData = dwData || ((enable && 1)<<15);
setmxcsr(dwData)
/*
	__asm__ __volatile__ (
		movzx	ecx, byte ptr enable
		and		ecx, 1
		shl		ecx, 15
		STMXCSR	dword ptr dwData
		mov		eax, dwData
		and		eax, ~(1<<15)	// clear FTZ bit
		or		eax, ecx	// set the FTZ bit
		mov		dwData, eax
		LDMXCSR	dword ptr dwData
		)*/
}




} //end System
} //end SGF
#endif //WIN32

