
#include "xmmintrin.h"
#include "sys/SMF_System.h"

#define MXCSR_DAZ (1 << 6)	/* Enable denormals are zero mode */
#define MXCSR_FTZ (1 << 15)	/* Enable flush to zero mode */

#define FXSAVE	(1 << 24)
#define SSE	(1 << 25)



namespace SMF{
namespace System{
static void __attribute__((constructor))
#ifndef __x86_64__
/* The i386 ABI only requires 4-byte stack alignment, so this is necessary
   to make sure the fxsave struct gets correct alignment.
   See PR27537 and PR28621.  */
__attribute__ ((force_align_arg_pointer))
#endif


CFPU::fpu_HasDazSupport (bool &result)
{
#ifndef __x86_64__
  /* All 64-bit targets have SSE and DAZ; only check them explicitly
     for 32-bit ones. */
  unsigned int eax, ebx, ecx, edx;
  #if defined(MASM_INTEL)
  /* See if we can use cpuid.  */
    __asm__ ("pushf{l|d}\n\t"
	   "pushf{l|d}\n\t"
	   "pop{l}\t%0\n\t"
	   "mov{l}\t{%0, %1|%1, %0}\n\t"
	   "xor{l}\t{%2, %0|%0, %2}\n\t"
	   "push{l}\t%0\n\t"
	   "popf{l|d}\n\t"
	   "pushf{l|d}\n\t"
	   "pop{l}\t%0\n\t"
	   "popf{l|d}\n\t"
	   : "=&r" (eax), "=&r" (ebx)
	   : "i" (0x00200000));
  #else
  /* See if we can use cpuid.  */
  asm volatile ("pushfl; pushfl; popl %0; movl %0,%1; xorl %2,%0;"
		"pushl %0; popfl; pushfl; popl %0; popfl"
		: "=&r" (eax), "=&r" (ebx)
		: "i" (0x00200000));
  #endif
  if (((eax ^ ebx) & 0x00200000) == 0)
    return;
 #if defined(MASM_INTEL)
   /* Check the highest input value for eax.  */
  asm volatile ("xchg  %1,ebx; cpuid; xchg  %1, ebx"
		: "=a" (eax), "=r" (ebx), "=c" (ecx), "=d" (edx)
		: "0" (0));

 #else
   /* Check the highest input value for eax.  */
  asm volatile ("xchgl %%ebx, %1; cpuid; xchgl %%ebx, %1"
		: "=a" (eax), "=r" (ebx), "=c" (ecx), "=d" (edx)
		: "0" (0));
#endif

  if (eax == 0)
    return;
 #if defined(MASM_INTEL)
 asm volatile ("xchg %1,ebx ; cpuid; xchg  %1,ebx"
		: "=a" (eax), "=r" (ebx), "=c" (ecx), "=d" (edx)
		: "0" (1));
#else
  asm volatile ("xchgl %%ebx, %1; cpuid; xchgl %%ebx, %1"
		: "=a" (eax), "=r" (ebx), "=c" (ecx), "=d" (edx)
		: "0" (1));
#endif
  if (edx & SSE)
    {
     // unsigned int mxcsr = __builtin_ia32_stmxcsr ();

      //mxcsr |= MXCSR_FTZ;

      if (edx & FXSAVE)
	{
	  /* Check if DAZ is available.  */
	  struct
	    {
	      unsigned short int cwd;
	      unsigned short int swd;
	      unsigned short int twd;
	      unsigned short int fop;
	      long int fip;
	      long int fcs;
	      long int foo;
	      long int fos;
	      long int mxcsr;
	      long int mxcsr_mask;
	      long int st_space[32];
	      long int xmm_space[32];
	      long int padding[56];
	    } __attribute__ ((aligned (16))) fxsave;

	  __builtin_memset (&fxsave, 0, sizeof (fxsave));

	  asm volatile ("fxsave %0" : "=m" (fxsave) : "m" (fxsave));

	  if (fxsave.mxcsr_mask & MXCSR_DAZ){
	    result =true;
	    //mxcsr |= MXCSR_DAZ;
	    }else {result =false;}
	}

      //__builtin_ia32_ldmxcsr (mxcsr);
    }
#else
  unsigned int mxcsr = __builtin_ia32_stmxcsr ();
  mxcsr |= MXCSR_DAZ | MXCSR_FTZ;
  __builtin_ia32_ldmxcsr (mxcsr);
#endif
}

static void __attribute__((constructor))
#ifndef __x86_64__
/* The i386 ABI only requires 4-byte stack alignment, so this is necessary
   to make sure the fxsave struct gets correct alignment.
   See PR27537 and PR28621.  */
__attribute__ ((force_align_arg_pointer))
#endif


CFPU::fpu_getFXSave (void *buffer)
{
#ifndef __x86_64__
  /* All 64-bit targets have SSE and DAZ; only check them explicitly
     for 32-bit ones. */
  unsigned int eax, ebx, ecx, edx;

#if defined(MASM_INTEL)
    __asm__ ("pushf{l|d}\n\t"
	   "pushf{l|d}\n\t"
	   "pop{l}\t%0\n\t"
	   "mov{l}\t{%0, %1|%1, %0}\n\t"
	   "xor{l}\t{%2, %0|%0, %2}\n\t"
	   "push{l}\t%0\n\t"
	   "popf{l|d}\n\t"
	   "pushf{l|d}\n\t"
	   "pop{l}\t%0\n\t"
	   "popf{l|d}\n\t"
	   : "=&r" (eax), "=&r" (ebx)
	   : "i" (0x00200000));
  #else
  /* See if we can use cpuid.  */
  asm volatile ("pushfl; pushfl; popl %0; movl %0,%1; xorl %2,%0;"
		"pushl %0; popfl; pushfl; popl %0; popfl"
		: "=&r" (eax), "=&r" (ebx)
		: "i" (0x00200000));
  #endif
  if (((eax ^ ebx) & 0x00200000) == 0)
    return;
#if defined(MASM_INTEL)
 asm volatile ("xchg %1,ebx ; cpuid; xchg  %1,ebx"
		: "=a" (eax), "=r" (ebx), "=c" (ecx), "=d" (edx)
		: "0" (1));
#else
  /* Check the highest input value for eax.  */
  asm volatile ("xchgl %%ebx, %1; cpuid; xchgl %%ebx, %1"
		: "=a" (eax), "=r" (ebx), "=c" (ecx), "=d" (edx)
		: "0" (0));
#endif
  if (eax == 0)
    return;
#if defined(MASM_INTEL)
  asm volatile ("xchg  %1,ebx; cpuid; xchg  %1,ebx"
		: "=a" (eax), "=r" (ebx), "=c" (ecx), "=d" (edx)
		: "0" (1));
#else
  asm volatile ("xchgl %%ebx, %1; cpuid; xchgl %%ebx, %1"
		: "=a" (eax), "=r" (ebx), "=c" (ecx), "=d" (edx)
		: "0" (1));
#endif
  if (edx & SSE)
    {
     // unsigned int mxcsr = __builtin_ia32_stmxcsr ();

      //mxcsr |= MXCSR_FTZ;

      if (edx & FXSAVE)
    {
	  /* Check if DAZ is available.  */
	  struct
	    {
	      unsigned short int cwd;
	      unsigned short int swd;
	      unsigned short int twd;
	      unsigned short int fop;
	      long int fip;
	      long int fcs;
	      long int foo;
	      long int fos;
	      long int mxcsr;
	      long int mxcsr_mask;
	      long int st_space[32];
	      long int xmm_space[32];
	      long int padding[56];
	    } __attribute__ ((aligned (16))) fxsave;

	  __builtin_memset (&fxsave, 0, sizeof (fxsave));

	  asm volatile ("fxsave %0" : "=m" (fxsave) : "m" (fxsave));
      __builtin_memcpy(buffer, &fxsave, sizeof(buffer));

	}

      //__builtin_ia32_ldmxcsr (mxcsr);
    }
#else
  unsigned int mxcsr = __builtin_ia32_stmxcsr ();
  mxcsr |= MXCSR_DAZ | MXCSR_FTZ;
  __builtin_ia32_ldmxcsr (mxcsr);
#endif
}

} //end Global
} //end SGF

