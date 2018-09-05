
#if defined(_WIN32) && !defined(_MSC_VER)
#include "sys/SGF_System.h"
#include "util/SGF_Heap.h"
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <linux/unistd.h>
#include <sys/syscall.h>
#include <sys/types.h>
#include <cpuid.h>


namespace SGF{
namespace System{
/*
==============================================================

	Clock ticks

==============================================================
*/

/*
================
getClockTicks
================
*/
#if defined(__i386__)
#if 0
__inline__ double getClockTicks(void)
{
    unsigned long long int x;
    __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
    return x;
}
#endif
__inline__ double getClockTicks(void)
{
   uint64_t x;
  __asm__ volatile ("rdtsc" : "=A" (x));
  return x;

}
#elif defined(__x86_64__)
#if 0
__inline__ double getClockTicks(void)
{
    unsigned hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}
#endif
__inline__ double getClockTicks(void)
{
  uint64_t a, d;
  __asm__ volatile ("rdtsc" : "=a" (a), "=d" (d));
  return (d<<32) | a;
}

#endif
#if 0
double getClockTicks( void ) {


	LARGE_INTEGER li;

	QueryPerformanceCounter( &li );
	return = (double ) li.LowPart + (double) 0xFFFFFFFF * li.HighPart;

//#else

	unsigned long lo, hi;
	uint64_t tick1; //unsigned 64 bit quantity
	asm volatile("pushl %%ebx \n\t");
	asm volatile("xor %eax, %eax\n\t"
		     "cpuid\n\t");

	asm volatile("rdtsc" : "=a" (lo), "=d" (hi));  //assembly code running the instruction
        asm volatile("popl %%ebx\n\t");                                                          //rdtsc

	tick1 = (((uint64_t)lo) | (((uint64_t)hi) << 32)); // calculating the tick value.
	return tick1; //(double ) lo + (double) 0xFFFFFFFF * hi;





}
#endif



/*
================
getClockTicksPerSecond
================
*/
double getClockTicksPerSecond( void ) {
	static double ticks = 0;
#if 0

	if ( !ticks ) {
		LARGE_INTEGER li;
		QueryPerformanceFrequency( &li );
		ticks = li.QuadPart;
	}



	if ( !ticks ) {
		HKEY hKey;
		LPBYTE ProcSpeed;
		DWORD buflen, ret;

		if ( !RegOpenKeyEx( HKEY_LOCAL_MACHINE, "HARDWARE\\DESCRIPTION\\System\\CentralProcessor\\0", 0, KEY_READ, &hKey ) ) {
			ProcSpeed = 0;
			buflen = sizeof( ProcSpeed );
			ret = RegQueryValueEx( hKey, "~MHz", NULL, NULL, (LPBYTE) &ProcSpeed, &buflen );
			// If we don't succeed, try some other spellings.
			if ( ret != ERROR_SUCCESS ) {
				ret = RegQueryValueEx( hKey, "~Mhz", NULL, NULL, (LPBYTE) &ProcSpeed, &buflen );
			}
			if ( ret != ERROR_SUCCESS ) {
				ret = RegQueryValueEx( hKey, "~mhz", NULL, NULL, (LPBYTE) &ProcSpeed, &buflen );
			}
			RegCloseKey( hKey );
			if ( ret == ERROR_SUCCESS ) {
				ticks = (double) ((unsigned long)ProcSpeed) * 1000000;
			}
		}
	}

#endif
	return ticks;
}


/*
==============================================================

	CPU

==============================================================
*/

/*
================
HasCPUID
================
*/


static bool HasCPUID( void ) {
unsigned int __eax, __ebx, __ecx, __edx;
#ifndef __x86_64__
  /* See if we can use cpuid.  On AMD64 we always can.  */
#if __GNUC__ >= 3
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
	   : "=&r" (__eax), "=&r" (__ebx)
	   : "i" (0x00200000));
#else
/* Host GCCs older than 3.0 weren't supporting Intel asm syntax
   nor alternatives in i386 code.  */
  __asm__ ("pushfl\n\t"
	   "pushfl\n\t"
	   "popl\t%0\n\t"
	   "movl\t%0, %1\n\t"
	   "xorl\t%2, %0\n\t"
	   "pushl\t%0\n\t"
	   "popfl\n\t"
	   "pushfl\n\t"
	   "popl\t%0\n\t"
	   "popfl\n\t"
	   : "=&r" (__eax), "=&r" (__ebx)
	   : "i" (0x00200000));
#endif

  if (!((__eax ^ __ebx) & 0x00200000))
    return false;
#endif

return true;
}




#define HT_BIT             0x10000000     // EDX[28]  Bit 28 is set if HT is supported
#define FAMILY_ID          0x0F00         // EAX[11:8] Bit 8-11 contains family processor ID.
#define PENTIUM4_ID        0x0F00
#define EXT_FAMILY_ID      0x0F00000      // EAX[23:20] Bit 20-23 contains extended family processor ID
#define NUM_LOGICAL_BITS   0x00FF0000     // EBX[23:16] Bit 16-23 in ebx contains the number of logical
                                          // processors per physical processor when execute cpuid with
                                          // eax set to 1

#define INITIAL_APIC_ID_BITS  0xFF000000  // EBX[31:24] Bits 24-31 (8 bits) return the 8-bit unique
                                          // initial APIC ID for the processor this code is running on.
                                          // Default value = 0xff if HT is not supported



#define _REG_EAX		0
#define _REG_EBX		1
#define _REG_ECX		2
#define _REG_EDX		3

/*
================
CPUID
================
*/
//
// func=0: Get vendor ID
// func=1: Processor Info and Feature Bits
// func=2: Cache and TLB Descriptor information
// func=3: Processor Serial Number
// func=80000000h: Get Highest Extended Function Supported
// func=80000001h: Extended Processor Info and Feature Bits
// func=80000002h,80000003h,80000004h: Processor Brand String
// func=80000005h: L1 Cache and TLB Identifiers
// func=80000006h: Extended L2 Cache Features
// func=80000007h: Advanced Power Management Information
// func=80000008h: Virtual and Physical address Sizes
/**
 * Calls cpuid with op and store results of eax,ebx,ecx,edx
 * \param func cpuid function (eax input)
 * \param regs[4]  content of eax, ebx,ecx,edx  after the call to cpuid
 */
static void CPUID( uint32_t func,  uint32_t regs[4])  //unsigned int 32 bits
{
  #if 0
  uint32_t regEAX;
  uint32_t regEBX;
  uint32_t regECX;
  uint32_t regEDX;


  // GCC won't allow us to clobber EBX since its used to store the GOT. So we need to
  // lie to GCC and backup/restore EBX without declaring it as clobbered.
#if defined(__i386__)
  asm volatile( "pushl %%ebx   \n\t"
#elif defined(__x86_64__)
    asm volatile( "push %%rbx   \n\t"
#endif

                "cpuid         \n\t"
                "movl %%ebx, %1\n\t"
#if defined(__i386__)
		"popl %%ebx    \n\t"
#elif defined(__x86_64__)
		"pop %%rbx    \n\t"
#endif

                : "=a"(regEAX), "=r"(regEBX), "=c"(regECX), "=d"(regEDX)
                : "a"(func)
                : "cc" );
  	regs[_REG_EAX] = regEAX;
	regs[_REG_EBX] = regEBX;
	regs[_REG_ECX] = regECX;
	regs[_REG_EDX] = regEDX;
#endif
__get_cpuid (func,
	     &regs[_REG_EAX], &regs[_REG_EBX],
	     &regs[_REG_ECX], &regs[_REG_EDX]);
}

/**
 * Retrieve the maximum function callable using cpuid
 */
uint32_t cpuid_maxcall()
{
  uint32_t regs[4];
  CPUID( 0, regs );
  return regs[_REG_EAX];
}

/**
 * Reference:
 * http://datasheets.chipdb.org/Intel/x86/CPUID/24161821.pdf
 * http://www.flounder.com/cpuid_explorer2.htm
 */
enum CpuidFeatures
{
  FPU   = 1<< 0, // Floating-Point Unit on-chip
  VME   = 1<< 1, // Virtual Mode Extension
  DE    = 1<< 2, // Debugging Extension
  PSE   = 1<< 3, // Page Size Extension
  TSC   = 1<< 4, // Time Stamp Counter
  MSR   = 1<< 5, // Model Specific Registers
  PAE   = 1<< 6, // Physical Address Extension
  MCE   = 1<< 7, // Machine Check Exception
  CX8   = 1<< 8, // CMPXCHG8 Instruction
  APIC  = 1<< 9, // On-chip APIC hardware
  SEP   = 1<<11, // Fast System Call
  MTRR  = 1<<12, // Memory type Range Registers
  PGE   = 1<<13, // Page Global Enable
  MCA   = 1<<14, // Machine Check Architecture
  CMOV  = 1<<15, // Conditional MOVe Instruction
  PAT   = 1<<16, // Page Attribute Table
  PSE36 = 1<<17, // 36bit Page Size Extension
  PSN   = 1<<18, // Processor Serial Number
  CLFSH = 1<<19, // CFLUSH Instruction
  DS    = 1<<21, // Debug Store
  ACPI  = 1<<22, // Thermal Monitor & Software Controlled Clock
  MMX   = 1<<23, // MultiMedia eXtension
  FXSR  = 1<<24, // Fast Floating Point Save & Restore
  SSE   = 1<<25, // Streaming SIMD Extension 1
  SSE2  = 1<<26, // Streaming SIMD Extension 2
  SS    = 1<<27, // Self Snoop
  HTT   = 1<<28, // Hyper Threading Technology
  TM    = 1<<29, // Thermal Monitor
  PBE   = 1<<31, // Pend Break Enabled
};
/**
 * This will retrieve the CPU features available
 * \return The content of the edx register containing available features
 */
uint32_t cpuid_features()
{
  uint32_t regs[4];
  CPUID( 1, regs );
  return regs[_REG_EDX];
}

/**
 * Reference:
 * http://datasheets.chipdb.org/Intel/x86/CPUID/24161821.pdf
 * http://www.flounder.com/cpuid_explorer2.htm
 */
enum CpuidExtendedFeatures
{
  SSE3  = 1<< 0, // Streaming SIMD Extension 3
  MW    = 1<< 3, // Mwait instruction
  CPL   = 1<< 4, // CPL-qualified Debug Store
  VMX   = 1<< 5, // VMX
  EST   = 1<< 7, // Enhanced Speed Test
  TM2   = 1<< 8, // Thermal Monitor 2
  L1    = 1<<10, // L1 Context ID
  CAE   = 1<<13, // CompareAndExchange 16B
};
/**
 * This will retrieve the extended CPU features available
 * \return The content of the ecx register containing available extended features
 */
uint32_t cpuid_extended_features()
{
  unsigned regs[4];
  CPUID( 1, regs);
  return regs[_REG_ECX];
}

/**
 * Retrieve the processor name.
 * \param name Preallocated string containing at least room for 13 characters. Will
 *             contain the name of the processor.
 */
void cpuid_procname( char* name )
{
  char pstring[16];
  char processorString[13];
  // get name of processor
  CPUID( 0, ( unsigned int * ) pstring );
	name[0] = pstring[4];
	name[1] = pstring[5];
	name[2] = pstring[6];
	name[3] = pstring[7];
	name[4] = pstring[12];
	name[5] = pstring[13];
	name[6] = pstring[14];
	name[7] = pstring[15];
	name[8] = pstring[8];
	name[9] = pstring[9];
	name[10] = pstring[10];
	name[11] = pstring[11];
	name[12] = 0;

}

/*
================
IsAMD
================
*/
static bool IsAMD( void ) {
	char pstring[16];
	char processorString[13]; //initialize it to zeros
	cpuid_procname( processorString );

	if ( strcmp( processorString, "AuthenticAMD" ) == 0 ) {
		return true;
	}
	return false;
}

/*
================
HasCMOV
================
*/
static bool HasCMOV( void ) {
	unsigned regs[4];

	// get CPU feature bits
	CPUID( 1, regs );

	// bit 15 of EDX denotes CMOV existence
	if ( regs[_REG_EDX] & ( 1 << 15 ) ) {
		return true;
	}
	return false;
}

/*
================
Has3DNow
================
*/
static bool Has3DNow( void ) {
	unsigned regs[4];

	// check AMD-specific functions
	CPUID( 0x80000000, regs );
	if ( regs[_REG_EAX] < 0x80000000 ) {
		return false;
	}

	// bit 31 of EDX denotes 3DNow! support
	CPUID( 0x80000001, regs );
	if ( regs[_REG_EDX] & ( 1 << 31 ) ) {
		return true;
	}

	return false;
}

/*
================
HasMMX
================
*/
static bool HasMMX( void ) {
	unsigned regs[4];

	// get CPU feature bits
	CPUID( 1, regs );

	// bit 23 of EDX denotes MMX existence
	if ( regs[_REG_EDX] & ( 1 << 23 ) ) {
		return true;
	}
	return false;
}

/*
================
HasSSE
================
*/
static bool HasSSE( void ) {
	unsigned regs[4];

	// get CPU feature bits
	CPUID( 1, regs );

	// bit 25 of EDX denotes SSE existence
	if ( regs[_REG_EDX] & ( 1 << 25 ) ) {
		return true;
	}
	return false;
}

/*
================
HasSSE2
================
*/
static bool HasSSE2( void ) {
	unsigned regs[4];

	// get CPU feature bits
	CPUID( 1, regs );

	// bit 26 of EDX denotes SSE2 existence
	if ( regs[_REG_EDX] & ( 1 << 26 ) ) {
		return true;
	}
	return false;
}

/*
================
HasSSE3
================
*/
static bool HasSSE3( void ) {
	unsigned regs[4];

	// get CPU feature bits
	CPUID( 1, regs );

	// bit 0 of ECX denotes SSE3 existence
	if ( regs[_REG_ECX] & ( 1 << 0 ) ) {
		return true;
	}
	return false;
}
/*
* Function "HTSupported"-detects whether Hyper-Threading technology is supported.
* This function will check to see if the CPU is a genuine Intel CPU and if bit 28 of the register EDX is set to 1.
* Note that the above conditions only indicate that Hyper-threading is supported, not enabled.
*/
unsigned int HTSupported(void)
{


	unsigned int Regedx= 0, Regeax= 0, VendorId[3] = {0, 0, 0};
	unsigned regs[4];
	// Get vendor id string
	CPUID( 0, regs );
	VendorId[0]=regs[_REG_EBX];
	VendorId[1]=regs[_REG_EDX];
	VendorId[2]=regs[_REG_ECX];
	// call cpuid with eax = 1
	CPUID( 0, regs );
	Regeax=regs[_REG_EAX]; // eax contains family processor type
	Regedx=regs[_REG_EDX]; // edx has info about the availability of hyper-Threading

	if (((Regeax & FAMILY_ID) ==  PENTIUM4_ID) ||
		(Regeax & EXT_FAMILY_ID))
	  if (VendorId[0] == 'uneG')
		if (VendorId[1] == 'Ieni')
			if (VendorId[2] == 'letn')
			  return(Regedx & HT_BIT);    // Genuine Intel with hyper-Threading technology

	return 0;    // Not genuine Intel processor

}

/*
================
LogicalProcPerPhysicalProc
================
*/
#define NUM_LOGICAL_BITS   0x00FF0000     // EBX[23:16] Bit 16-23 in ebx contains the number of logical
                                          // processors per physical processor when execute cpuid with
                                          // eax set to 1

/*Function "LogicalProcPerPhysicalProc"-returns the number of logical processors associating
with one physica l processor.

*/
static unsigned char LogicalProcPerPhysicalProc( void ) {
	unsigned int regebx = 0;
	unsigned regs[4];
	CPUID( 1, regs );
	regebx = regs[_REG_EBX];
	return (unsigned char) ((regebx & NUM_LOGICAL_BITS) >> 16);
}

/*
================
GetAPIC_ID
================
*/
#define INITIAL_APIC_ID_BITS  0xFF000000  // EBX[31:24] Bits 24-31 (8 bits) return the 8-bit unique
                                          // initial APIC ID for the processor this code is running on.
                                          // Default value = 0xff if HT is not supported
/*Function "GetAPIC_ID" to return the ID of the logical processor. To identify whether Hyper-Threading
technology is enabled or not is to check if there exist more than one logical processor ID per physical processor ID.
*/
static unsigned char GetAPIC_ID( void ) {
	unsigned int regebx = 0;
	unsigned regs[4];
	CPUID( 1, regs );
	regebx = regs[_REG_EBX];
	return (unsigned char) ((regebx & INITIAL_APIC_ID_BITS) >> 24);
}

static pid_t gettid( void )
{
    return syscall( __NR_gettid );
}

/*
================
CPUCount

	logicalNum is the number of logical CPU per physical CPU
    physicalNum is the total number of physical processor
	returns one of the HT_* flags
	http://software.intel.com/en-us/articles/counting-physical-and-logical-32-bit-processors
================
*/
#define HT_NOT_CAPABLE			0
#define HT_ENABLED			1
#define HT_DISABLED			2
#define HT_SUPPORTED_NOT_ENABLED	3
#define HT_CANNOT_DETECT		4


int CPUCount( int &logicalNum, int &physicalNum ) {
	int statusFlag;

	physicalNum = -1;
	logicalNum = 1;
	statusFlag = HT_NOT_CAPABLE;


	// Number of physical processors in a non-Intel system
	// or in a 32-bit Intel system with Hyper-Threading technology disabled
	physicalNum = sysconf(_SC_NPROCESSORS_ONLN);
	if (physicalNum < 1)
	{
	fprintf(stderr, "Could not determine number of CPUs online:\n%s\n",
	strerror (errno));

	}
	unsigned char HT_Enabled = 0;
if (HTSupported())
	{
	logicalNum = LogicalProcPerPhysicalProc();
	// Will Verify HT flags
	/*
	 * Basically, this is done using a bitmask. This means that there is an integer (say, 32 bits),
	 * and if the first bit == 1, then that process is allowed to run on processor 1. If the second bit == 1,
	 * then that process is allowed to run on processor 2. etc.
	 * So by default, the affinity bitmask = 1...111 (32 times.) This means that "The process may run on processor 1, 2, 3, ..., 32."
	 * Of course, if you only have 2 cores, then the extra 30 bits won't apply.
	 * However, if you set that bitmask to be: 0...010 then only "processor 2" is allowed to execute that process.
	 * http://linux.die.net/man/2/sched_setaffinity
	 */

	if ( logicalNum >= 1 ) {	// > 1 doesn't mean HT is enabled in the BIOS

		cpu_set_t dwProcessAffinity;        /* Define your cpu_set bit mask. */
		cpu_set_t dwSystemAffinity;
		cpu_set_t dwAffinityMask;
		CPU_ZERO(&dwProcessAffinity);       /* Initialize it all to 0, i.e. no CPUs selected. */
		CPU_ZERO(&dwSystemAffinity);
		CPU_ZERO(&dwAffinityMask);
		// Calculate the appropriate  shifts and mask based on the
		// number of logical processors.

		unsigned char i = 1, PHY_ID_MASK  = 0xFF, PHY_ID_SHIFT = 0;

		while( i < logicalNum ) {
			i *= 2;
 			PHY_ID_MASK  <<= 1;
			PHY_ID_SHIFT++;
		}

		size_t  mask_len;

		// GET CURRENT PROCESS AFFINITY MASK
		int ret = sched_getaffinity(0, mask_len, &dwProcessAffinity);
		if (ret) {
		  printf("sched_get_affinity returned %d, exiting.\n", ret);
		  statusFlag = HT_CANNOT_DETECT;
		  physicalNum = -1;
		  return statusFlag;
		}
		printf("current process's affinity: %l bytes mask.\n", mask_len);
		printf("CPU_COUNT() of set:    %d\n", CPU_COUNT(&dwProcessAffinity));

		int count =1;
		CPU_SET(count, &dwAffinityMask);

		while ( count != 0 && count <=CPU_COUNT(&dwProcessAffinity)  ) {
			// Check if this CPU is available
			if ( CPU_ISSET(count, &dwProcessAffinity)) {
			  if (sched_setaffinity( gettid(), sizeof( cpu_set_t ), &dwAffinityMask ))
			  {
				perror( "sched_setaffinity" );
				return NULL;
			  }else
				{
					unsigned char APIC_ID, LOG_ID, PHY_ID;

					Util::rest( 1 ); // Give OS time to switch CPU

					APIC_ID = GetAPIC_ID();
					LOG_ID  = APIC_ID & ~PHY_ID_MASK;
					PHY_ID  = APIC_ID >> PHY_ID_SHIFT;

					if ( LOG_ID != 0 ) {
						HT_Enabled = 1;
					}
				}
			}
			CPU_SET(count++, &dwAffinityMask);
		}

		// Reset the processor affinity
		sched_setaffinity(0, sizeof(cpu_set_t), &dwProcessAffinity); /* Set affinity of tihs process to */
                                                  /* the defined mask, i.e. only 7. */
		if ( logicalNum == 1 ) {  // Normal P4 : HT is disabled in hardware
			statusFlag = HT_DISABLED;
		} else {
			if ( HT_Enabled ) {
				// Total physical processors in a Hyper-Threading enabled system.
				physicalNum /= logicalNum;
				statusFlag = HT_ENABLED;
			} else {
				statusFlag = HT_SUPPORTED_NOT_ENABLED;
			}
		}
	}
  }
  else
  {
	// Processors do not have Hyper-Threading technology
	statusFlag = HT_NOT_CAPABLE;
        logicalNum = 1;
  }
  return statusFlag;
}


/*
================
HasHTT
================
*/
static bool HasHTT( void ) {
	unsigned regs[4];
	int logicalNum, physicalNum, HTStatusFlag;

	// get CPU feature bits
	CPUID( 1, regs );

	// bit 28 of EDX denotes HTT existence
	if ( !( regs[_REG_EDX] & ( 1 << 28 ) ) ) {
		return false;
	}

	HTStatusFlag = CPUCount( logicalNum, physicalNum );
	if ( HTStatusFlag != HT_ENABLED ) {
		return false;
	}
	return true;
}

/*
================
HasHTT
================
*/
static bool HasDAZ( void ) {
	unsigned char FXSaveArea[512] __attribute__ ((aligned (16)));
	unsigned char *FXArea = FXSaveArea;
	SGF_Dword dwMask = 0;
	unsigned regs[4];

	// get CPU feature bits
	CPUID( 1, regs );

	// bit 24 of EDX denotes support for FXSAVE
	if ( !( regs[_REG_EDX] & ( 1 << 24 ) ) ) {
		return false;
	}

	memset( FXArea, 0, sizeof( FXSaveArea ) );
	asm volatile("fxsave %0" : : "m" (FXArea));

	dwMask = *(SGF_Dword *)&FXArea[28];	// Read the MXCSR Mask
	return ( ( dwMask & ( 1 << 6 ) ) == ( 1 << 6 ) );	// Return if the DAZ bit is set
}
/*
================
Sys_GetCPUId
================
*/
static cpuid_t getCPUId( void ) {
	int flags;

	// verify we're at least a Pentium or 486 with CPUID support
	if ( !HasCPUID() ) {
		return CPUID_UNSUPPORTED;
	}

	// check for an AMD
	if ( IsAMD() ) {
		flags = CPUID_AMD;
	} else {
		flags = CPUID_INTEL;
	}

	// check for Multi Media Extensions
	if ( HasMMX() ) {
		flags |= CPUID_MMX;
	}

	// check for 3DNow!
	if ( Has3DNow() ) {
		flags |= CPUID_3DNOW;
	}

	// check for Streaming SIMD Extensions
	if ( HasSSE() ) {
		flags |= CPUID_SSE | CPUID_FTZ;
	}

	// check for Streaming SIMD Extensions 2
	if ( HasSSE2() ) {
		flags |= CPUID_SSE2;
	}

	// check for Streaming SIMD Extensions 3 aka Prescott's New Instructions
	if ( HasSSE3() ) {
		flags |= CPUID_SSE3;
	}

	// check for Hyper-Threading Technology
	if ( HasHTT() ) {
		flags |= CPUID_HTT;
	}

	// check for Conditional Move (CMOV) and fast floating point comparison (FCOMI) instructions
	if ( HasCMOV() ) {
		flags |= CPUID_CMOV;
	}

	// check for Denormals-Are-Zero mode
	if ( HasDAZ() ) {
		flags |= CPUID_DAZ;
	}

	return (cpuid_t)flags;
}





/*
================
Sys_GetProcessorString
================
*/
const char *getProcessorString( void ) {
	  cpuid_t			cpuid;
	  string string;
		//
	// CPU type
	//
	if (true /* !CMyString::Icmp( win32.sys_cpustring.GetString(), "detect" )*/ ) {


		Debug::debug(Debug::error,__FUNCTION__)<< getClockTicksPerSecond() / 1000000.0f << " MHz "<<endl;

		cpuid = getCPUId();

		string.clear();

		if ( cpuid & CPUID_AMD ) {
			string += "AMD CPU";
		} else if ( cpuid & CPUID_INTEL ) {
			string += "Intel CPU";
		} else if ( cpuid & CPUID_UNSUPPORTED ) {
			string += "unsupported CPU";
		} else {
			string += "generic CPU";
		}

		string += " with ";
		if ( cpuid & CPUID_MMX ) {
			string += "MMX & ";
		}
		if ( cpuid & CPUID_3DNOW ) {
			string += "3DNow! & ";
		}
		if ( cpuid & CPUID_SSE ) {
			string += "SSE & ";
		}
		if ( cpuid & CPUID_SSE2 ) {
            string += "SSE2 & ";
		}
		if ( cpuid & CPUID_SSE3 ) {
			string += "SSE3 & ";
		}
		if ( cpuid & CPUID_HTT ) {
			string += "HTT & ";
		}
//Todo:  StripTrailing
	//s	string.StripTrailing(string.c_str(), " & " );
	//s	string.StripTrailing( " with " );

	} /*else {
		Debug::debug(Debug::math,__FUNCTION__) << "forcing CPU type to " );
		idLexer src( win32.sys_cpustring.GetString(), CMyString::Length( win32.sys_cpustring.GetString() ), "sys_cpustring" );
		idToken token;

		int id = CPUID_NONE;
		while( src.ReadToken( &token ) ) {
			if ( token.Icmp( "generic" ) == 0 ) {
				id |= CPUID_GENERIC;
			} else if ( token.Icmp( "intel" ) == 0 ) {
				id |= CPUID_INTEL;
			} else if ( token.Icmp( "amd" ) == 0 ) {
				id |= CPUID_AMD;
			} else if ( token.Icmp( "mmx" ) == 0 ) {
				id |= CPUID_MMX;
			} else if ( token.Icmp( "3dnow" ) == 0 ) {
				id |= CPUID_3DNOW;
			} else if ( token.Icmp( "sse" ) == 0 ) {
				id |= CPUID_SSE;
			} else if ( token.Icmp( "sse2" ) == 0 ) {
				id |= CPUID_SSE2;
			} else if ( token.Icmp( "sse3" ) == 0 ) {
				id |= CPUID_SSE3;
			} else if ( token.Icmp( "htt" ) == 0 ) {
				id |= CPUID_HTT;
			}
		}
		if ( id == CPUID_NONE ) {
			Debug::debug(Debug::error,__FUNCTION__)<< "WARNING: unknown sys_cpustring: ", win32.sys_cpustring.GetString() );
			id = CPUID_GENERIC;
		}
		win32.cpuid = (cpuid_t) id;
	}
	*/
	return string.c_str();
}



/*
================
Sys_GetProcessorId
================
*/
cpuid_t getProcessorId( void ) {
    return getCPUId();
}
#if 0
////////////////////////////////////////////////////////////////////////////////////////
//  DLL and CLIPBOARD
//////////////////////////////
/*
================
Sys_GetClipboardData
================
*/
char *Sys_GetClipboardData( void ) {
	char *data = NULL;
	char *cliptext;

	if ( OpenClipboard( NULL ) != 0 ) {
		HANDLE hClipboardData;

		if ( ( hClipboardData = GetClipboardData( CF_TEXT ) ) != 0 ) {
			if ( ( cliptext = (char *)GlobalLock( hClipboardData ) ) != 0 ) {
				data = (char *)Mem_Alloc( GlobalSize( hClipboardData ) + 1 );
				strcpy( data, cliptext );
				GlobalUnlock( hClipboardData );

				strtok( data, "\n\r\b" );
			}
		}
		CloseClipboard();
	}
	return data;
}

/*
================
Sys_SetClipboardData
================
*/
void Sys_SetClipboardData( const char *string ) {
	HGLOBAL HMem;
	char *PMem;

	// allocate memory block
	HMem = (char *)::GlobalAlloc( GMEM_MOVEABLE | GMEM_DDESHARE, strlen( string ) + 1 );
	if ( HMem == NULL ) {
		return;
	}
	// lock allocated memory and obtain a pointer
	PMem = (char *)::GlobalLock( HMem );
	if ( PMem == NULL ) {
		return;
	}
	// copy text into allocated memory block
	lstrcpy( PMem, string );
	// unlock allocated memory
	::GlobalUnlock( HMem );
	// open Clipboard
	if ( !OpenClipboard( 0 ) ) {
		::GlobalFree( HMem );
		return;
	}
	// remove current Clipboard contents
	EmptyClipboard();
	// supply the memory handle to the Clipboard
	SetClipboardData( CF_TEXT, HMem );
	HMem = 0;
	// close Clipboard
	CloseClipboard();
}

/*
========================================================================

DLL Loading

========================================================================
*/
static const int		MAX_OSPATH					= 256;
#define INTSIGNBITNOTSET(i)		((~((const unsigned long)(i))) >> 31)
/*
================
CMyString::IcmpPath
================
*/
static int IcmpPath( const char *s1, const char *s2 ) {
	int c1, c2, d;

#if 0
//#if !defined( _WIN32 )
	Debug::debug(Debug::math,__FUNCTION__) << "WARNING: IcmpPath used on a case-sensitive filesystem?\n" );
#endif

	do {
		c1 = *s1++;
		c2 = *s2++;

		d = c1 - c2;
		while( d ) {
			if ( c1 <= 'Z' && c1 >= 'A' ) {
				d += ('a' - 'A');
				if ( !d ) {
					break;
				}
			}
			if ( c1 == '\\' ) {
				d += ('/' - '\\');
				if ( !d ) {
					break;
				}
			}
			if ( c2 <= 'Z' && c2 >= 'A' ) {
				d -= ('a' - 'A');
				if ( !d ) {
					break;
				}
			}
			if ( c2 == '\\' ) {
				d -= ('/' - '\\');
				if ( !d ) {
					break;
				}
			}
			// make sure folders come first
			while( c1 ) {
				if ( c1 == '/' || c1 == '\\' ) {
					break;
				}
				c1 = *s1++;
			}
			while( c2 ) {
				if ( c2 == '/' || c2 == '\\' ) {
					break;
				}
				c2 = *s2++;
			}
			if ( c1 && !c2 ) {
				return -1;
			} else if ( !c1 && c2 ) {
				return 1;
			}
			// same folder depth so use the regular compare
			return ( INTSIGNBITNOTSET( d ) << 1 ) - 1;
		}
	} while( c1 );

	return 0;
}

/*
=====================
Sys_DLL_Load
=====================
*/
int Sys_DLL_Load( const char *dllName ) {
	HINSTANCE	libHandle;
	libHandle = LoadLibrary( dllName );
	if ( libHandle ) {
		// since we can't have LoadLibrary load only from the specified path, check it did the right thing
		char loadedPath[ MAX_OSPATH ];
		GetModuleFileName( libHandle, loadedPath, sizeof( loadedPath ) - 1 );
		if ( IcmpPath( dllName, loadedPath ) ) {
			Debug::debug(Debug::error,__FUNCTION__)<< "ERROR: LoadLibrary "<< dllName<<" wants to load :"<< loadedPath <<endl;
			Sys_DLL_Unload( (int)libHandle );
			return 0;
		}
	}
	return (int)libHandle;
}

/*
=====================
Sys_DLL_GetProcAddress
=====================
*/
void *Sys_DLL_GetProcAddress( int dllHandle, const char *procName ) {
	return GetProcAddress( (HINSTANCE)dllHandle, procName );
}

/*
=====================
Sys_DLL_Unload
=====================
*/
void Sys_DLL_Unload( int dllHandle ) {
	if ( !dllHandle ) {
		return;
	}
	if ( FreeLibrary( (HINSTANCE)dllHandle ) == 0 ) {
		int lastError = GetLastError();
		LPVOID lpMsgBuf;
		FormatMessage(
			FORMAT_MESSAGE_ALLOCATE_BUFFER,
		    NULL,
			lastError,
			MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), // Default language
			(LPTSTR) &lpMsgBuf,
			0,
			NULL
		);
		//s Sys_Error( "Sys_DLL_Unload: FreeLibrary failed - %s (%d)", lpMsgBuf, lastError );
	}
}
#endif //end if0
} //end namespace System
} //end SGF
#if 0  //example
#include <cstdlib>
#include <iostream>

int main( int, char** )
{
  char procname[13];
  cpuid_procname(procname);
  std::cout << "Processor name: " << procname << std::endl;
  std::cout << std::endl;
  std::cout << "Max cpuid call: " << cpuid_maxcall() << std::endl;
  std::cout << std::endl;
  std::cout << "Processor features:" << std::endl;
  std::cout << "  FPU   = " << std::boolalpha << (bool)(cpuid_features() & FPU  ) << std::endl;
  std::cout << "  VME   = " << std::boolalpha << (bool)(cpuid_features() & VME  ) << std::endl;
  std::cout << "  DE    = " << std::boolalpha << (bool)(cpuid_features() & DE   ) << std::endl;
  std::cout << "  PSE   = " << std::boolalpha << (bool)(cpuid_features() & PSE  ) << std::endl;
  std::cout << "  TSC   = " << std::boolalpha << (bool)(cpuid_features() & TSC  ) << std::endl;
  std::cout << "  MSR   = " << std::boolalpha << (bool)(cpuid_features() & MSR  ) << std::endl;
  std::cout << "  PAE   = " << std::boolalpha << (bool)(cpuid_features() & PAE  ) << std::endl;
  std::cout << "  MCE   = " << std::boolalpha << (bool)(cpuid_features() & MCE  ) << std::endl;
  std::cout << "  CX8   = " << std::boolalpha << (bool)(cpuid_features() & CX8  ) << std::endl;
  std::cout << "  APIC  = " << std::boolalpha << (bool)(cpuid_features() & APIC ) << std::endl;
  std::cout << "  SEP   = " << std::boolalpha << (bool)(cpuid_features() & SEP  ) << std::endl;
  std::cout << "  MTRR  = " << std::boolalpha << (bool)(cpuid_features() & MTRR ) << std::endl;
  std::cout << "  PGE   = " << std::boolalpha << (bool)(cpuid_features() & PGE  ) << std::endl;
  std::cout << "  MCA   = " << std::boolalpha << (bool)(cpuid_features() & MCA  ) << std::endl;
  std::cout << "  CMOV  = " << std::boolalpha << (bool)(cpuid_features() & CMOV ) << std::endl;
  std::cout << "  PAT   = " << std::boolalpha << (bool)(cpuid_features() & PAT  ) << std::endl;
  std::cout << "  PSE36 = " << std::boolalpha << (bool)(cpuid_features() & PSE36) << std::endl;
  std::cout << "  PSN   = " << std::boolalpha << (bool)(cpuid_features() & PSN  ) << std::endl;
  std::cout << "  CLFSH = " << std::boolalpha << (bool)(cpuid_features() & CLFSH) << std::endl;
  std::cout << "  DS    = " << std::boolalpha << (bool)(cpuid_features() & DS   ) << std::endl;
  std::cout << "  ACPI  = " << std::boolalpha << (bool)(cpuid_features() & ACPI ) << std::endl;
  std::cout << "  MMX   = " << std::boolalpha << (bool)(cpuid_features() & MMX  ) << std::endl;
  std::cout << "  FXSR  = " << std::boolalpha << (bool)(cpuid_features() & FXSR ) << std::endl;
  std::cout << "  SSE   = " << std::boolalpha << (bool)(cpuid_features() & SSE  ) << std::endl;
  std::cout << "  SSE2  = " << std::boolalpha << (bool)(cpuid_features() & SSE2 ) << std::endl;
  std::cout << "  SS    = " << std::boolalpha << (bool)(cpuid_features() & SS   ) << std::endl;
  std::cout << "  HTT   = " << std::boolalpha << (bool)(cpuid_features() & HTT  ) << std::endl;
  std::cout << "  TM    = " << std::boolalpha << (bool)(cpuid_features() & TM   ) << std::endl;
  std::cout << std::endl;
  std::cout << "Processor extended features:" << cpuid_extended_features() << std::endl;
  std::cout << "  SSE3 = " << std::boolalpha << (bool)(cpuid_extended_features() & SSE3) << std::endl;
  std::cout << "  MW   = " << std::boolalpha << (bool)(cpuid_extended_features() & MW  ) << std::endl;
  std::cout << "  CPL  = " << std::boolalpha << (bool)(cpuid_extended_features() & CPL ) << std::endl;
  std::cout << "  VMX  = " << std::boolalpha << (bool)(cpuid_extended_features() & VMX ) << std::endl;
  std::cout << "  EST  = " << std::boolalpha << (bool)(cpuid_extended_features() & EST ) << std::endl;
  std::cout << "  TM2  = " << std::boolalpha << (bool)(cpuid_extended_features() & TM2 ) << std::endl;
  std::cout << "  L1   = " << std::boolalpha << (bool)(cpuid_extended_features() & L1  ) << std::endl;
  std::cout << "  CAE  = " << std::boolalpha << (bool)(cpuid_extended_features() & CAE ) << std::endl;

  return EXIT_SUCCESS;
}


#endif
#endif //_WIN32
