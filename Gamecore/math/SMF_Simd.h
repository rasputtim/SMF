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

#ifndef _SMF_MATH_SIMD_H__
#define _SMF_MATH_SIMD_H__
#include "../SMF_Config.h"
#include "../util/SMF_UtilStructs.h"
#include "../util/SMF_ArgsCmdLine.h"
#include "../math/SMF_Math.h"



namespace SMF{
using namespace Util;


namespace GEO{
class CPlane;
class CVertex;
class CAABBox;
}
namespace MATH{
class CQuaternion;
class CVec3D;
class CJointQuaternion;
class CMatJoint3x4;
extern unsigned start_ClockCount;
extern unsigned end_ClockCount;
extern unsigned saved_ebx;

#if defined(_MSC_VER)
/*example: how to get correct timestamp information:
cpuid ; force all previous instructions to complete
rdtsc ; read time stamp counter
mov time, eax ; move counter into variable
fdiv ; floating-point divide (medindo quantos ciclos demora para fazer a divisão)
cpuid ; wait for FDIV to complete before RDTSC
rdtsc ; read time stamp counter
sub eax, time ; find the difference

*/
//TODO: testar com getClockTicks()
#define RecordTime(var) \
__asm cpuid \
__asm rdtsc \
__asm mov var, eax

#define StartRecordTime \
__asm mov SMF::MATH::saved_ebx, ebx \
RecordTime(SMF::MATH::start_ClockCount)

#define StopRecordTime \
RecordTime(SMF::MATH::end_ClockCount) \
__asm mov ebx, SMF::MATH::saved_ebx
//TODO testar com System::startRecordingTimeStamp(start)
#define StartRecordTimeLocal( start )			\
	__asm mov SMF::MATH::saved_ebx, ebx				\
	__asm xor eax, eax						\
	RecordTime(start)						\
	__asm xor eax, eax						\
	__asm cpuid
//\todo: testar com stopRecordingTimStamp(double &end)
#define StopRecordTimeLocal( end )				\
	__asm xor eax, eax		 				\
	RecordTime(end)							\
	__asm mov ebx, SMF::MATH::saved_ebx				\
	__asm xor eax, eax						\
	__asm cpuid

#define GetBest( start, end, best )			\
	if ( !best || (end - start < best) ) {	\
		best = end - start;					\
	}
#endif

#if defined(__GNUC__) //for GNU
#define STRINGIZE(X) #X
#if defined(MASM_INTEL)

#define RecordTime(var) \
"cpuid\n" \
"rdtsc\n" \
"mov  "STRINGIZE(var)",eax \n"

#define StartRecordTime \
asm("mov %0, ebx\n"      \
RecordTime(%1) \
::"m"(SMF::MATH::saved_ebx),"m"(SMF::MATH::start_ClockCount):);


#define StopRecordTime \
asm(            \
    RecordTime(%1) \
"mov  %0,ebx\n"::"m"(SMF::MATH::saved_ebx),"m"(SMF::MATH::end_ClockCount):);




#define StartRecordTimeLocal( start )    \
asm( "mov  %0, ebx\n"				\
	"xor eax, eax\n"						\
	RecordTime(%1)						\
	"xor eax, eax\n"						\
	"cpuid\n"                       \
	::"m"(SMF::MATH::saved_ebx),"m"(start) :);


#define StopRecordTimeLocal( end )      \
asm( 	\
	"xor eax, eax\n"						\
	RecordTime(%1)						\
	"mov  ebx,%0\n"              \
	"xor eax, eax\n"						\
	"cpuid\n"                       \
	::"m"(SMF::MATH::saved_ebx),"m"(end) :);  \


#else

#define RecordTime(var) \
"cpuid\n" \
"rdtsc\n" \
"movl %%eax,"STRINGIZE(var)"\n"

#define StartRecordTime \
asm("movl %%ebx, %0\n"      \
RecordTime(%1) \
::"m"(SMF::MATH::saved_ebx),"m"(SMF::MATH::start_ClockCount):);


#define StopRecordTime \
asm(            \
    RecordTime(%1) \
"movl %%ebx, %0\n"::"m"(SMF::MATH::saved_ebx),"m"(SMF::MATH::end_ClockCount):);




#define StartRecordTimeLocal( start )    \
asm( "movl %%ebx, %0\n"				\
	"xor %%eax, %%eax\n"						\
	RecordTime(%1)						\
	"xor %%eax, %%eax\n"						\
	"cpuid\n"                       \
	::"m"(SMF::MATH::saved_ebx),"m"(start) :);


#define StopRecordTimeLocal( end )  \
	asm( 	\
	"xor %%eax, %%eax\n"						\
	RecordTime(%1)						\
	"movl %0, %%ebx\n"              \
	"xor %%eax, %%eax\n"						\
	"cpuid\n"                       \
	::"m"(SMF::MATH::saved_ebx),"m"(end) :);  \


#endif
#endif // defined  __GNUC__


// Computes the inverse of a 4x4 matrix via direct cofactor expansion.
/// Returns the determinant of the original matrix, and zero on failure.
#define MAT_COFACTOR(mat, i, j) \
	_mm_sub_ps(_mm_mul_ps(_mm_shuffle_ps(mat[2], mat[1], _MM_SHUFFLE(j,j,j,j)), \
	           shuffle1_ps(_mm_shuffle_ps(mat[3], mat[2], _MM_SHUFFLE(i,i,i,i)), _MM_SHUFFLE(2,0,0,0))), \
	           _mm_mul_ps(shuffle1_ps(_mm_shuffle_ps(mat[3], mat[2], _MM_SHUFFLE(j,j,j,j)), _MM_SHUFFLE(2,0,0,0)), \
	           _mm_shuffle_ps(mat[2], mat[1], _MM_SHUFFLE(i,i,i,i))))
#define _mm_transpose_matrix_intel(row0, row1, row2, row3) \
	sf_m128 tmp0, tmp1, tmp2, tmp3; \
	tmp0 = _mm_unpacklo_ps(row0, row1); \
	tmp2 = _mm_unpacklo_ps(row2, row3); \
	tmp1 = _mm_unpackhi_ps(row0, row1); \
	tmp3 = _mm_unpackhi_ps(row2, row3); \
	row0 = _mm_movelh_ps(tmp0, tmp2); \
	row1 = _mm_movehl_ps(tmp2, tmp0); \
	row2 = _mm_movelh_ps(tmp1, tmp3); \
	row3 = _mm_movehl_ps(tmp3, tmp1);
// Computes the inverse of a 4x4 matrix via direct cofactor expansion.
/// Returns the determinant of the original matrix, and zero on failure.
#define MAT_COFACTOR(mat, i, j) \
	_mm_sub_ps(_mm_mul_ps(_mm_shuffle_ps(mat[2], mat[1], _MM_SHUFFLE(j,j,j,j)), \
	           shuffle1_ps(_mm_shuffle_ps(mat[3], mat[2], _MM_SHUFFLE(i,i,i,i)), _MM_SHUFFLE(2,0,0,0))), \
	           _mm_mul_ps(shuffle1_ps(_mm_shuffle_ps(mat[3], mat[2], _MM_SHUFFLE(j,j,j,j)), _MM_SHUFFLE(2,0,0,0)), \
	           _mm_shuffle_ps(mat[2], mat[1], _MM_SHUFFLE(i,i,i,i))))

const sf_u32 andMaskOne = 0xFFFFFFFF;
const float andMaskOneF = MATH::ReinterpretAsFloat(andMaskOne);

/// A SSE mask register with x = y = z = 0xFFFFFFFF and w = 0x0.
extern const sf_m128 sseMaskXYZ;

extern const sf_m128 sseSignMask ; // -0.f = 1 << 31
extern const sf_m128 sseSignMask3; // -0.f = 1 << 31
extern const sf_m128 sseSignMask4;

extern const __m256 sseSignMask256; // -0.f = 1 << 31

///\todo Benchmark which one is better!
//#define negate_ps(x) _mm_xor_ps(x, sseSignMask)
#define negate_ps(x) _mm_sub_ps(_mm_setzero_ps(), x)
#define negate3_ps(x) _mm_xor_ps(x, sseSignMask3)

/**
\brief Returns the lowest element of the given sse register as a float.
\note When compiling with /arch:SSE or newer, it is expected that this function is a no-op "cast", since
the resulting float is represented in an XMM register as well. Check the disassembly to confirm!
**/
float M128_TO_FLOAT(sf_m128 sse);

/**
\brief Returns a SSE variable with the given float f in the lowest index. The three higher indices are set to zero.
\note When compiling with /arch:SSE or newer, it is expected that this function is a no-op "cast" if the given
float is already in a register, since it will lie in an XMM register already. Check the disassembly to confirm!
\note Detected on VS2010 32-bit + AVX that this generates a vmovss+vxorps+vmovss instruction triple!
\note Never use this function if you need to generate a 4-vector [f,f,f,f]. Instead, use _mm_set1_ps(&f), which
generates a vmovss+vhufps and no redundant vxorps+vmovss!
**/

sf_m128 FLOAT_TO_M128(float f);


sf_m128 SET_PS(float um,float dois, float tres, float quatro);

/**
\brief Given four scalar SS FP registers, packs the four values into a single SP FP register.
**/
SMF_INLINE sf_m128 pack_4ss_to_ps(sf_m128 x, sf_m128 y, sf_m128 z, const sf_m128 &w)
{
	sf_m128 xy = _mm_movelh_ps(x, y); // xy = [ _, y, _, x]
	sf_m128 zw = _mm_movelh_ps(z, w); // zw = [ _, w, _, z]
	return _mm_shuffle_ps(xy, zw, _MM_SHUFFLE(2, 0, 2, 0)); // ret = [w, z, y, x]
}
// If mask[i] == 0, then output index i from a, otherwise mask[i] must be 0xFFFFFFFF, and output index i from b.
SMF_INLINE  sf_m128 cmov_ps(sf_m128 a, sf_m128 b, sf_m128 mask)
{
//#ifdef MATH_SSE41 // SSE 4.1 offers conditional copying between registers with the blendvps instruction.
//	return _mm_blendv_ps(a, b, mask);
//#else // If not on SSE 4.1, use conditional masking.
	b = _mm_and_ps(mask, b); // Where mask is 1, output b.
	a = _mm_andnot_ps(mask, a); // Where mask is 0, output a.
	return _mm_or_ps(a, b);
//#endif
}

//FOR SSE2
#define set_ps_hex2(w, z, y, x) _mm_castsi128_ps(_mm_set_epi32(w, z, y, x))
#define set1_ps_hex2(x) _mm_castsi128_ps(_mm_set1_epi32(x))
//FOR SSE
#define set_ps_hex(w, z, y, x) _mm_set_ps(ReinterpretAsFloat(w), ReinterpretAsFloat(z), ReinterpretAsFloat(y), ReinterpretAsFloat(x))
#define set1_ps_hex(x) _mm_set1_ps(ReinterpretAsFloat(x))

//for SSE
#define abs_ps(x) _mm_andnot_ps(sseSignMask, x)
//for AVX
#define abs_ps256(x) _mm256_andnot_ps(sseSignMask256, x)


//===============================================================
/**
 * \class 	CProcClock

 *
 * \ingroup SMF_Math
 *
 * \if pt_br
 * \brief	Medidas de ciclos de clock do processador
 *
 * \elseif us_en
 * \brief 	Processor Timestamp clock cicles measurements
 *
 * \endif
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CProcClock {
public:
    static void printClocks( const char *string, int dataCount, int clocks, int otherClocks = 0 );
    static void getBaseClocks();
};

class CSIMDProcessor;

/**
 * \class CSIMD
 *
 * \ingroup SMF_Math
 *
 * \brief Implementa um processador SIMD
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *

 * \note http://softpixel.com/~cwright/programming/simd/cpuid.php

===============================================================================

	Single Instruction Multiple Data (SIMD)

	For optimal use data should be aligned on a 16 sf_u8 boundary.
	All CSIMDProcessor routines are thread safe.

===============================================================================
A sigla SIMD (Single Instruction, Multiple Data), descreve um método de operação de
computadores com várias unidades operacionais em computação paralela. Neste modo, a mesma
instrução é aplicada
simultaneamente a diversos dados para produzir mais resultados.
O modelo SIMD é adequado para o tratamento de conjuntos regulares de dados, como as matrizes
e vetores. Esse tipo de máquina opera aplicando uma única instrução a um conjunto de elementos
de um vetor. Sendo uma máquina que aplique a n elementos uma determinada instrução e o vetor t
contenha os elementos a serem processados, t terá todos seus elementos calculados n vezes mais
rápido que uma máquina SISD na mesma tarefa.
http://pt.wikipedia.org/wiki/SIMD
http://softpixel.com/~cwright/programming/simd/cpuid.php
*/
class SMF_API CSIMD {
public:
	/**
	\brief inicializa o processador SIMD genérico
	\note iniciliza o processador genérico
	**/
	static void			init();
	/**
	\brief verifica se o processador generico foi inicializado
	\return verdadeiro caso o processador genérico tenha sido inicializado, falso em caso contrário
		**/
	static bool			isInitialized();
	/**
	\brief inicializa o processador SIMD
	\note descobre qual o processador atual e constrói a classe equivalente
	**/
	static void			initProcessor( const char *module, bool forceGeneric );
	/**
	\brief verifica se o processador específico foi inicializado
	\return verdadeiro caso o processador específico tenha sido inicializado, falso em caso contrário
	**/
	static bool			isProcessorInitialized();
	/**
	\brief retorna o SIMD processor.
	\note se não houver sido inicializado, inicializa e retorna
	**/
	static  CSIMDProcessor * CSIMD::getProcessor();
	/**
	\brief retorna o Generic processor.
	\note se não houver sido inicializado, inicializa e retorna
	**/
	static  CSIMDProcessor * CSIMD::getGenProcessor();


	/**
	\brief wrapper to CMatXD
	**/
	static bool			globalUseSIMD();
	/**
	\brief wrapper to CMatXD
	**/
	static void         initHeap();
	/**
	\brief deleta as estruturas de dados utilizadas
	**/
	static void			shutdown();
	/**
	\brief programas de teste de performance
	**/
	static void			Test_f( const class CCMDLineArgs &args );

	/**
	\brief programas de teste de performance de CVec3D
	**/
	static void			Test3D_f( const class CCMDLineArgs &args );

	static void         TesteLenght3D(const class CCMDLineArgs &args);
    /**
	\brief programas de teste de performance de CVec4D
	**/
	static void			Test4D_f( const class CCMDLineArgs &args );
	/**
	\brief programas de teste de performance de CMat4D
	**/
	static void			TestMat4D_f();
	/**
	\brief programas de teste de performance de CQuaternion
	**/
	static void			TestQuat4D_f();
	/**
	\brief programas de teste de performance de CPlane
	**/
	static void			TestPlane4D_f();
	/**
	\brief programas de teste das rotinas de calculo trigonométrico
	**/
	static void         SSE_TestTrigonometry();

};


/*
===============================================================================

	virtual base class for different SIMD processors
	//http://www.popoloski.com/posts/sse_move_instructions/
===============================================================================
*/


using namespace GEO;
class CVec2D;
class CVec3D;
class CVec4D;
class CVec5D;
class CVec6D;
class CVecXD;
class CMat2D;
class CMat3D;
class CMat4D;
class CMat5D;
class CMat6D;
class CMatXD;

class CQuaternion;
//s struct dominantTri_s;

const int MIXBUFFER_SAMPLES = 4096;

typedef enum {
	SPEAKER_LEFT = 0,
	SPEAKER_RIGHT,
	SPEAKER_CENTER,
	SPEAKER_LFE,
	SPEAKER_BACKLEFT,
	SPEAKER_BACKRIGHT
} speakerLabel;


/**
 * \class CSIMDProcessor
 *
 * \ingroup SMF_Math
 *
 * \brief classe base virtual base paradiferentes processadores SIMD
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *

 * \note http://www.popoloski.com/posts/sse_move_instructions/
 **/
class SMF_API CSIMDProcessor {
public:
	CSIMDProcessor():
	cpuid(CPUID_NONE){  }
	CSIMDProcessor(cpuid_t id):
	cpuid(id){  }
	~CSIMDProcessor(){};
	virtual const char * VPCALL	getName() const = 0;
	/**
	\brief   dst[i] = constant + src[i];
	\param dst ponteiro para um buffer de floats onde será armazenado o resultado
	\param src ponteiro para um buffer de floats
	\param constant constant que será somada a src
	\param count qtde de floats em cada buffer - (i)
	**/
	virtual void VPCALL add( float *dst,			const float constant,	const float *src,		const int count ) = 0;

	/**
	\brief   dst[i] = src0[i] + src1[i];
	\param dst ponteiro para um buffer de floats onde será armazenado o resultado
	\param src0 ponteiro para um buffer de floats
	\param src1 ponteiro para um buffer de floats
	\param count qtde de floats em cada buffer - (i)
	**/
	virtual void VPCALL add( float *dst,			const float *src0,		const float *src1,		const int count ) = 0;
	/**
	\brief   dst[i] = constant - src[i];
	\param dst ponteiro para um buffer de floats onde será armazenado o resultado
	\param src ponteiro para um buffer de floats
	\param constant constant que será somada a src
	\param count qtde de floats em cada buffer - (i)
	**/
	virtual void VPCALL sub( float *dst,			const float constant,	const float *src,		const int count ) = 0;
	/**
	\brief   dst[i] = src0[i] - src1[i];
	\param dst ponteiro para um buffer de floats onde será armazenado o resultado
	\param src0 ponteiro para um buffer de floats
	\param src1 ponteiro para um buffer de floats
	\param count qtde de floats em cada buffer - (i)
	**/
	virtual void VPCALL sub( float *dst,			const float *src0,		const float *src1,		const int count ) = 0;
	/**
	\brief   dst[i] = constant * src[i];
	\param dst ponteiro para um buffer de floats onde será armazenado o resultado
	\param src ponteiro para um buffer de floats
	\param constant constant que será somada a src
	\param count qtde de floats em cada buffer - (i)
	**/
	virtual void VPCALL mul( float *dst,			const float constant,	const float *src,		const int count ) = 0;
	/**
	\brief   dst[i] = src0[i] * src1[i];
	\param dst ponteiro para um buffer de floats onde será armazenado o resultado
	\param src0 ponteiro para um buffer de floats
	\param src1 ponteiro para um buffer de floats
	\param count qtde de floats em cada buffer - (i)
	**/
	virtual void VPCALL mul( float *dst,			const float *src0,		const float *src1,		const int count ) = 0;
	/**
	\brief   dst[i] = constant / src[i];
	\param dst ponteiro para um buffer de floats onde será armazenado o resultado
	\param src ponteiro para um buffer de floats
	\param constant constant que será somada a src
	\param count qtde de floats em cada buffer - (i)
	**/
	virtual void VPCALL div( float *dst,			const float constant,	const float *src,		const int count ) = 0;
	/**
	\brief   dst[i] = src0[i] / src1[i];
	\param dst ponteiro para um buffer de floats onde será armazenado o resultado
	\param src0 ponteiro para um buffer de floats
	\param src1 ponteiro para um buffer de floats
	\param count qtde de floats em cada buffer - (i)
	**/
	virtual void VPCALL div( float *dst,			const float *src0,		const float *src1,		const int count ) = 0;

	/**
	\brief   dst[i] += constant * src[i];
	\param dst ponteiro para um buffer de floats onde será armazenado o resultado
	\param src ponteiro para um buffer de floats
	\param constant constant que será somada a src
	\param count qtde de floats em cada buffer - (i)
	**/
	virtual void VPCALL mulAdd( float *dst,			const float constant,	const float *src,		const int count ) = 0;
	/**
	\brief   dst[i] += constant * src[i];
	\param dst ponteiro para um buffer de floats onde será armazenado o resultado
	\param src0 ponteiro para um buffer de floats
	\param src1 ponteiro para um buffer de floats
	\param count qtde de floats em cada buffer - (i)
	**/
	virtual void VPCALL mulAdd( float *dst,			const float *src0,		const float *src1,		const int count ) = 0;
	/**
	\brief   dst[i] -= constant * src[i];
	\param dst ponteiro para um buffer de floats onde será armazenado o resultado
	\param src ponteiro para um buffer de floats
	\param constant constant que será somada a src
	\param count qtde de floats em cada buffer - (i)
	**/
	virtual void VPCALL mulSub( float *dst,			const float constant,	const float *src,		const int count ) = 0;
	/**
	\brief  dst[i] -= constant * src[i];
	\param dst ponteiro para um buffer de floats onde será armazenado o resultado
	\param src0 ponteiro para um buffer de floats
	\param src1 ponteiro para um buffer de floats
	\param count qtde de floats em cada buffer - (i)
	**/
	virtual void VPCALL mulSub( float *dst,			const float *src0,		const float *src1,		const int count ) = 0;
	/**
	\brief  produto de vetor por uma constante. dst[i] = constant * src[i](Vector);
	\param dst ponteiro para o buffer que armazenará os resultados
	\param constant constante que ser´multiplicada pelo vetor
	\param src buffer contendo vetoresn3D
	\param count número de oerações a realizar. (i)
	**/
	virtual	void VPCALL dot( float *dst,			const CVec3D &constant,	const CVec3D *src,		const int count ) = 0;
	virtual	void VPCALL dot( float *dst,			const CVec3D &constant,	const CPlane *src,		const int count ) = 0;
	virtual void VPCALL dot( float *dst,			const CVec3D &constant,	const CVertex *src,	const int count ) = 0;
	virtual	void VPCALL dot( float *dst,			const CPlane &constant,const CVec3D *src,		const int count ) = 0;
	virtual	void VPCALL dot( float *dst,			const CPlane &constant,const CPlane *src,		const int count ) = 0;
	virtual void VPCALL dot( float *dst,			const CPlane &constant,const CVertex *src,	const int count ) = 0;
	virtual	void VPCALL dot( float *dst,			const CVec3D *src0,		const CVec3D *src1,		const int count ) = 0;
	virtual void VPCALL dot( float &dot,			const float *src1,		const float *src2,		const int count ) = 0;

	virtual	void VPCALL cmpGT( sf_u8 *dst,			const float *src0,		const float constant,	const int count ) = 0;
	virtual	void VPCALL cmpGT( sf_u8 *dst,			const sf_u8 bitNum,		const float *src0,		const float constant,	const int count ) = 0;
	virtual	void VPCALL cmpGE( sf_u8 *dst,			const float *src0,		const float constant,	const int count ) = 0;
	virtual	void VPCALL cmpGE( sf_u8 *dst,			const sf_u8 bitNum,		const float *src0,		const float constant,	const int count ) = 0;
	virtual	void VPCALL cmpLT( sf_u8 *dst,			const float *src0,		const float constant,	const int count ) = 0;
	virtual	void VPCALL cmpLT( sf_u8 *dst,			const sf_u8 bitNum,		const float *src0,		const float constant,	const int count ) = 0;
	virtual	void VPCALL cmpLE( sf_u8 *dst,			const float *src0,		const float constant,	const int count ) = 0;
	virtual	void VPCALL cmpLE( sf_u8 *dst,			const sf_u8 bitNum,		const float *src0,		const float constant,	const int count ) = 0;

	virtual	void VPCALL minMax( float &min,			float &max,				const float *src,		const int count ) = 0;
	virtual	void VPCALL minMax( CVec2D &min,		CVec2D &max,			const CVec2D *src,		const int count ) = 0;
	virtual	void VPCALL minMax( CVec3D &min,		CVec3D &max,			const CVec3D *src,		const int count ) = 0;
	virtual	void VPCALL minMax( CVec3D &min,		CVec3D &max,			const CVertex *src,	const int count ) = 0;
	virtual	void VPCALL minMax( CVec3D &min,		CVec3D &max,			const CVertex *src,	const int *indexes,		const int count ) = 0;

	virtual	void VPCALL clamp( float *dst,			const float *src,		const float min,		const float max,		const int count ) = 0;
	virtual	void VPCALL clampMin( float *dst,		const float *src,		const float min,		const int count ) = 0;
	virtual	void VPCALL clampMax( float *dst,		const float *src,		const float max,		const int count ) = 0;

	virtual void VPCALL memCopy( void *dst,			const void *src,		const int count ) = 0;
	virtual void VPCALL memSet( void *dst,			const int val,			const int count ) = 0;

	// these assume 16 sf_u8 aligned and 16 sf_u8 padded memory
	/**
	\note Este método assume que haja um alinhamento e preenchimento de memória de 16 bytes
	\note O malloc do compilador não consegue garantir o alinhamento da memória,
	\note apresentando muitos erros nos procedimentos.
	\note deve ser utilizado o CHeap para garantir o alinhamento da memória.
	**/
	virtual void VPCALL zero16( float *dst,			const int count ) = 0;
	/**
	\note Este método assume que haja um alinhamento e preenchimento de memória de 16 bytes
	\note O malloc do compilador não consegue garantir o alinhamento da memória,
	\note apresentando muitos erros nos procedimentos.
	\note deve ser utilizado o CHeap para garantir o alinhamento da memória.
	**/
	virtual void VPCALL negate16( float *dst,		const int count ) = 0;
	/**
	\note Este método assume que haja um alinhamento e preenchimento de memória de 16 bytes
	\note O malloc do compilador não consegue garantir o alinhamento da memória,
	\note apresentando muitos erros nos procedimentos.
	\note deve ser utilizado o CHeap para garantir o alinhamento da memória.
	**/
	virtual void VPCALL copy16( float *dst,			const float *src,		const int count ) = 0;
	/**
	\note Este método assume que haja um alinhamento e preenchimento de memória de 16 bytes
	\note O malloc do compilador não consegue garantir o alinhamento da memória,
	\note apresentando muitos erros nos procedimentos.
	\note deve ser utilizado o CHeap para garantir o alinhamento da memória.
	**/
	virtual void VPCALL add16( float *dst,			const float *src1,		const float *src2,		const int count ) = 0;
	/**
	\note Este método assume que haja um alinhamento e preenchimento de memória de 16 bytes
	\note O malloc do compilador não consegue garantir o alinhamento da memória,
	\note apresentando muitos erros nos procedimentos.
	\note deve ser utilizado o CHeap para garantir o alinhamento da memória.
	**/
	virtual void VPCALL sub16( float *dst,			const float *src1,		const float *src2,		const int count ) = 0;
	/**
	\note Este método assume que haja um alinhamento e preenchimento de memória de 16 bytes
	\note O malloc do compilador não consegue garantir o alinhamento da memória,
	\note apresentando muitos erros nos procedimentos.
	\note deve ser utilizado o CHeap para garantir o alinhamento da memória.
	**/
	virtual void VPCALL mul16( float *dst,			const float *src1,		const float constant,	const int count ) = 0;
	/**
	\note Este método assume que haja um alinhamento e preenchimento de memória de 16 bytes
	\note O malloc do compilador não consegue garantir o alinhamento da memória,
	\note apresentando muitos erros nos procedimentos.
	\note deve ser utilizado o CHeap para garantir o alinhamento da memória.
	**/
	virtual void VPCALL addAssign16( float *dst,	const float *src,		const int count ) = 0;
	/**
	\note Este método assume que haja um alinhamento e preenchimento de memória de 16 bytes
	\note O malloc do compilador não consegue garantir o alinhamento da memória,
	\note apresentando muitos erros nos procedimentos.
	\note deve ser utilizado o CHeap para garantir o alinhamento da memória.
	**/
	virtual void VPCALL subAssign16( float *dst,	const float *src,		const int count ) = 0;
	/**
	\note Este método assume que haja um alinhamento e preenchimento de memória de 16 bytes
	\note O malloc do compilador não consegue garantir o alinhamento da memória,
	\note apresentando muitos erros nos procedimentos.
	\note deve ser utilizado o CHeap para garantir o alinhamento da memória.
	**/
	virtual void VPCALL mulAssign16( float *dst,	const float constant,	const int count ) = 0;

	// CMatXD operations

	/**
	\brief Inverte a matriz 4x4 passada
	\note Este método assume que haja um alinhamento e preenchimento de memória de 16 bytes
	\note O malloc do compilador não consegue garantir o alinhamento da memória,
	\note apresentando muitos erros nos procedimentos.
	\note deve ser utilizado o CHeap para garantir o alinhamento da memória.
	\return false se não houver inversa
	**/

	virtual bool VPCALL matX_inverse_4x4(float* src)=0;

	/**
	\note Este método assume que haja um alinhamento e preenchimento de memória de 16 bytes
	\note O malloc do compilador não consegue garantir o alinhamento da memória,
	\note apresentando muitos erros nos procedimentos.
	\note deve ser utilizado o CHeap para garantir o alinhamento da memória.
	**/
	virtual void VPCALL matX_MultiplyVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec ) = 0;
	/**
	\note Este método assume que haja um alinhamento e preenchimento de memória de 16 bytes
	\note O malloc do compilador não consegue garantir o alinhamento da memória,
	\note apresentando muitos erros nos procedimentos.
	\note deve ser utilizado o CHeap para garantir o alinhamento da memória.
	**/
	virtual void VPCALL matX_MultiplyAddVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec ) = 0;
	/**
	\note Este método assume que haja um alinhamento e preenchimento de memória de 16 bytes
	\note O malloc do compilador não consegue garantir o alinhamento da memória,
	\note apresentando muitos erros nos procedimentos.
	\note deve ser utilizado o CHeap para garantir o alinhamento da memória.
	**/
	virtual void VPCALL matX_MultiplySubVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec ) = 0;
	/**
	\note Este método assume que haja um alinhamento e preenchimento de memória de 16 bytes
	\note O malloc do compilador não consegue garantir o alinhamento da memória,
	\note apresentando muitos erros nos procedimentos.
	\note deve ser utilizado o CHeap para garantir o alinhamento da memória.
	**/
	virtual void VPCALL matX_TransposeMultiplyVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec ) = 0;
	/**
	\note Este método assume que haja um alinhamento e preenchimento de memória de 16 bytes
	\note O malloc do compilador não consegue garantir o alinhamento da memória,
	\note apresentando muitos erros nos procedimentos.
	\note deve ser utilizado o CHeap para garantir o alinhamento da memória.
	**/
	virtual void VPCALL matX_TransposeMultiplyAddVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec ) = 0;
	/**
	\note Este método assume que haja um alinhamento e preenchimento de memória de 16 bytes
	\note O malloc do compilador não consegue garantir o alinhamento da memória,
	\note apresentando muitos erros nos procedimentos.
	\note deve ser utilizado o CHeap para garantir o alinhamento da memória.
	**/
	virtual void VPCALL matX_TransposeMultiplySubVecX( CVecXD &dst, const CMatXD &mat, const CVecXD &vec ) = 0;
	/**
	\note Este método assume que haja um alinhamento e preenchimento de memória de 16 bytes
	\note O malloc do compilador não consegue garantir o alinhamento da memória,
	\note apresentando muitos erros nos procedimentos.
	\note deve ser utilizado o CHeap para garantir o alinhamento da memória.
	**/
	virtual void VPCALL matX_MultiplyMatX( CMatXD &dst, const CMatXD &m1, const CMatXD &m2 ) = 0;
	/**
	\note Este método assume que haja um alinhamento e preenchimento de memória de 16 bytes
	\note O malloc do compilador não consegue garantir o alinhamento da memória,
	\note apresentando muitos erros nos procedimentos.
	\note deve ser utilizado o CHeap para garantir o alinhamento da memória.
	**/
	virtual void VPCALL matX_TransposeMultiplyMatX( CMatXD &dst, const CMatXD &m1, const CMatXD &m2 ) = 0;
	/**
	\note Este método assume que haja um alinhamento e preenchimento de memória de 16 bytes
	\note O malloc do compilador não consegue garantir o alinhamento da memória,
	\note apresentando muitos erros nos procedimentos.
	\note deve ser utilizado o CHeap para garantir o alinhamento da memória.
	**/
	virtual void VPCALL matX_LowerTriangularsolve( const CMatXD &L, float *x, const float *b, const int n, int skip = 0 ) = 0;
	/**
	\note Este método assume que haja um alinhamento e preenchimento de memória de 16 bytes
	\note O malloc do compilador não consegue garantir o alinhamento da memória,
	\note apresentando muitos erros nos procedimentos.
	\note deve ser utilizado o CHeap para garantir o alinhamento da memória.
	**/
	virtual void VPCALL matX_LowerTriangularsolveTranspose( const CMatXD &L, float *x, const float *b, const int n ) = 0;
	/**
	\note Este método assume que haja um alinhamento e preenchimento de memória de 16 bytes
	\note O malloc do compilador não consegue garantir o alinhamento da memória,
	\note apresentando muitos erros nos procedimentos.
	\note deve ser utilizado o CHeap para garantir o alinhamento da memória.
	**/
	virtual bool VPCALL matX_LDLTFactor( CMatXD &mat, CVecXD &invDiag, const int n ) = 0;

	// rendering
	virtual void VPCALL blendJoints( CJointQuaternion *joints, const CJointQuaternion *blendJoints, const float lerp, const int *index, const int numJoints ) = 0;
	virtual void VPCALL convertJointQuatsToJointMats( CMatJoint3x4 *jointMats, const CJointQuaternion *jointQuats, const int numJoints ) = 0;
	virtual void VPCALL convertJointMatsToJointQuats( CJointQuaternion *jointQuats, const CMatJoint3x4 *jointMats, const int numJoints ) = 0;
	virtual void VPCALL transformJoints( CMatJoint3x4 *jointMats, const int *parents, const int firstJoint, const int lastJoint ) = 0;
	virtual void VPCALL untransformJoints( CMatJoint3x4 *jointMats, const int *parents, const int firstJoint, const int lastJoint ) = 0;
	virtual void VPCALL transformVerts( CVertex *verts, const int numVerts, const CMatJoint3x4 *joints, const CVec4D *weights, const int *index, const int numWeights ) = 0;
	virtual void VPCALL tracePointCull( sf_u8 *cullBits, sf_u8 &totalOr, const float radius, const CPlane *planes, const CVertex *verts, const int numVerts ) = 0;
	virtual void VPCALL decalPointCull( sf_u8 *cullBits, const CPlane *planes, const CVertex *verts, const int numVerts ) = 0;
	virtual void VPCALL overlayPointCull( sf_u8 *cullBits, CVec2D *texCoords, const CPlane *planes, const CVertex *verts, const int numVerts ) = 0;
	virtual void VPCALL deriveTriPlanes( CPlane *planes, const CVertex *verts, const int numVerts, const int *indexes, const int numIndexes ) = 0;
	virtual void VPCALL  deriveTangents( CPlane *planes, CVertex *verts, const int numVerts, const int *indexes, const int numIndexes ) = 0;
	//virtual void VPCALL deriveUnsmoothedTangents( CVertex *verts, const dominantTri_s *dominantTris, const int numVerts ) = 0;
	virtual void VPCALL  normalizeTangents( CVertex *verts, const int numVerts ) = 0;
	virtual void VPCALL  createTextureSpaceLightVectors( CVec3D *lightVectors, const CVec3D &lightOrigin, const CVertex *verts, const int numVerts, const int *indexes, const int numIndexes ) = 0;
	virtual void VPCALL  createSpecularTextureCoords( CVec4D *texCoords, const CVec3D &lightOrigin, const CVec3D &viewOrigin, const CVertex *verts, const int numVerts, const int *indexes, const int numIndexes ) = 0;
	virtual int  VPCALL  createShadowCache( CVec4D *vertexCache, int *vertRemap, const CVec3D &lightOrigin, const CVertex *verts, const int numVerts ) = 0;
	virtual int  VPCALL  createVertexProgramShadowCache( CVec4D *vertexCache, const CVertex *verts, const int numVerts ) = 0;

	// sound mixing
	virtual void VPCALL  upSamplePCMTo44kHz( float *dest, const short *pcm, const int numSamples, const int kHz, const int numChannels ) = 0;
	virtual void VPCALL  upSampleOGGTo44kHz( float *dest, const float * const *ogg, const int numSamples, const int kHz, const int numChannels ) = 0;
	virtual void VPCALL  mixSoundTwoSpeakerMono( float *mixBuffer, const float *samples, const int numSamples, const float lastV[2], const float currentV[2] ) = 0;
	virtual void VPCALL  mixSoundTwoSpeakerStereo( float *mixBuffer, const float *samples, const int numSamples, const float lastV[2], const float currentV[2] ) = 0;
	virtual void VPCALL  mixSoundSixSpeakerMono( float *mixBuffer, const float *samples, const int numSamples, const float lastV[6], const float currentV[6] ) = 0;
	virtual void VPCALL  mixSoundSixSpeakerStereo( float *mixBuffer, const float *samples, const int numSamples, const float lastV[6], const float currentV[6] ) = 0;
	virtual void VPCALL  mixedSoundToSamples( short *samples, const float *mixBuffer, const int numSamples ) = 0;

	//=================================


virtual void VPCALL vector3D_Sum(CVec3D* pOut, const CVec3D* pIn) = 0;
virtual void VPCALL vector3D_SumOf(CVec3D* pOut, const CVec3D* pIn1, const CVec3D* pIn2) = 0;
virtual void VPCALL vector3D_Diff(CVec3D* pLeft, CVec3D* pRight) = 0;
virtual void VPCALL vector3D_DiffOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight) = 0;
virtual void VPCALL vector3D_Scale(CVec3D* pOut, float scalar) = 0;
virtual void VPCALL vector3D_ScaleOf(CVec3D* pOut, const CVec3D* pIn, float scalar) = 0;
virtual float VPCALL vector3D_Dot(const CVec3D* pSrc1, const CVec3D* pSrc2) = 0;
virtual float VPCALL vector3D_LengthSq(const CVec3D* pVec) = 0;
virtual float VPCALL vector3D_Length(const CVec3D* pVec) = 0;
virtual void VPCALL vector3D_Cross(CVec3D* pLeft, const CVec3D* pRight) = 0;
virtual void VPCALL vector3D_CrossOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight) = 0;
virtual void VPCALL vector3D_Normalize(CVec3D* pVec) = 0;
virtual void VPCALL vector3D_NormalizeOf(CVec3D* pOut, const CVec3D* pVec) = 0;
virtual float VPCALL vector3D_Distance(const CVec3D* pVec1, const CVec3D* pVec2) = 0;

virtual void VPCALL vector3D_AlignedSum(CVec3D* pOut, const CVec3D* pIn) = 0;
virtual void VPCALL vector3D_AlignedSumOf(CVec3D* pOut, const CVec3D* pIn1, const CVec3D* pIn2) = 0;
virtual void VPCALL vector3D_AlignedDiff(CVec3D* pLeft, CVec3D* pRight) = 0;
virtual void VPCALL vector3D_AlignedDiffOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight) = 0;
virtual void VPCALL vector3D_AlignedScale(CVec3D* pOut, float scalar) = 0;
virtual void VPCALL vector3D_AlignedScaleOf(CVec3D* pOut, const CVec3D* pIn, float scalar) = 0;
virtual float VPCALL vector3D_AlignedDot(const CVec3D* pSrc1, const CVec3D* pSrc2) = 0;
virtual float VPCALL vector3D_AlignedLengthSq(const CVec3D* pVec) = 0;
virtual float VPCALL vector3D_AlignedLength(const CVec3D* pVec) = 0;
virtual void VPCALL vector3D_AlignedCross(CVec3D* pLeft, const CVec3D* pRight) = 0;
virtual void VPCALL vector3D_AlignedCrossOf(CVec3D* pOut, const CVec3D* pLeft, const CVec3D* pRight) = 0;
virtual void VPCALL vector3D_AlignedNormalize(CVec3D* pVec) = 0;
virtual void VPCALL vector3D_AlignedNormalizeOf(CVec3D* pOut, const CVec3D* pVec) = 0;
virtual float VPCALL vector3D_AlignedDistance(const CVec3D* pVec1, const CVec3D* pVec2) = 0;

virtual float VPCALL vector4D_Dot(const CVec4D* pSrc4D1, const CVec4D* pSrc4D2) = 0;
virtual void VPCALL vector4D_Sum(CVec4D* pOut, const CVec4D* pIn) = 0;
virtual void VPCALL vector4D_SumOf(CVec4D* pOut, const CVec4D* pIn1, const CVec4D* pIn2) = 0;
virtual void VPCALL vector4D_Diff(CVec4D* pLeft, CVec4D* pRight) = 0;
virtual void VPCALL vector4D_DiffOf(CVec4D* pOut, const CVec4D* pLeft, const CVec4D* pRight) = 0;
virtual void VPCALL vector4D_Scale(CVec4D* pOut, float scalar) = 0;
virtual void VPCALL vector4D_ScaleOf(CVec4D* pOut, const CVec4D* pIn, float scalar) = 0;
virtual float VPCALL vector4D_LengthSq(const CVec4D* pVec) = 0;
virtual float VPCALL vector4D_Length(const CVec4D* pVec) = 0;
virtual void VPCALL vector4D_Normalize(CVec4D* pVec) = 0;
virtual void VPCALL vector4D_NormalizeOf(CVec4D* pOut, const CVec4D* pVec) = 0;
virtual float VPCALL vector4D_Distance(const CVec4D* pVec1, const CVec4D* pVec2) = 0;

	virtual float VPCALL vector4D_AlignedDot(const CVec4D* pSrc4D1, const CVec4D* pSrc4D2) = 0;
	virtual void VPCALL vector4D_AlignedSum(CVec4D* pOut, const CVec4D* pIn) = 0;
	virtual void VPCALL vector4D_AlignedSumOf(CVec4D* pOut, const CVec4D* pIn1, const CVec4D* pIn2) = 0;
	virtual void VPCALL vector4D_AlignedDiff(CVec4D* pLeft, CVec4D* pRight) = 0;
	virtual void VPCALL vector4D_AlignedDiffOf(CVec4D* pOut, const CVec4D* pLeft, const CVec4D* pRight) = 0;
	virtual void VPCALL vector4D_AlignedScale(CVec4D* pOut, float scalar) = 0;
	virtual void VPCALL vector4D_AlignedScaleOf(CVec4D* pOut, const CVec4D* pIn, float scalar) = 0;
	virtual float VPCALL vector4D_AlignedLengthSq(const CVec4D* pVec) = 0;
	virtual float VPCALL vector4D_AlignedLength(const CVec4D* pVec) = 0;
	virtual void VPCALL vector4D_AlignedNormalize(CVec4D* pVec) = 0;
	virtual void VPCALL vector4D_AlignedNormalizeOf(CVec4D* pOut, const CVec4D* pVec) = 0;
	virtual float VPCALL vector4D_AlignedDistance(const CVec4D* pVec1, const CVec4D* pVec2) = 0;
//Trigonometry=====================

	virtual float VPCALL  invSqrt( float x ) = 0;
	virtual void  VPCALL  InvSqrt4( float x[4] ) = 0;
	virtual float VPCALL  sinZeroHalfPI( float x ) = 0;
	virtual void  VPCALL  sin4ZeroHalfPI( float entrada[4], float saida[4] ) = 0;
	virtual float VPCALL  sin( float a ) = 0;
	virtual void  VPCALL  sin4( float entrada[4], float saida[4] ) = 0;
	virtual float VPCALL  cosZeroHalfPI( float x ) = 0;
	virtual void  VPCALL  cos4ZeroHalfPI( float entrada[4], float saida[4] ) = 0;
	virtual float VPCALL  cos( float a ) = 0;
	virtual void  VPCALL  cos4( float entrada[4], float saida[4] ) = 0;
	virtual void  VPCALL  sincos( float a, float &sin, float &cos ) = 0;
	virtual void  VPCALL  sincos4( float a[4], float sin[4], float cos[4] ) = 0;
	virtual float VPCALL  aTanPositive( float y, float x ) = 0;
	virtual void  VPCALL  aTan4Positive( float y[4], float x[4], float resultado[4] ) = 0;
	virtual float VPCALL  atan( float y, float x ) = 0;
	virtual void  VPCALL  aTan4( float y[4], float x[4], float resultado[4] ) = 0;
	//===================

//======CMatriz==========================================



/**
\brief soma pMat e pIn e guarda o resultado em pMat
\param pMat Matriz que será somada a pIn e guardará o resultado da soma
\param pIn Matriz que será somada a pMat
\note Sums pMat and pIn and stores the result in pMat.
**/
virtual void  VPCALL mat4D_Sum(CMat4D* pMat, const CMat4D* pIn)=0;
/**
\brief soma pIn1 e pIn2 e armazena o resultado em pMat
\brief Sums pIn1 and pIn2 and stores the result in pMat.
\param pMat Matriz que guardará o resultado da soma
\param pIn1 Matriz que será somada a pIn2
\param pIn2 Matriz que será somada a pIn1
**/
virtual void  VPCALL mat4D_SumOf(CMat4D* pMat, const CMat4D* pIn1, const CMat4D* pIn2)=0;
/**
\brief subtrai pIn de pMat e armazena o resultado em pMat. (pMat = pMat - pIn)
\brief Subtracts pIn from pMat and stores the result in pMat.
\param pMat matriz que será utilizada para subtrair pIn e armazenará o resultado da subtração
\param pIn matriz que será subtrída de pMat
**/
virtual void  VPCALL mat4D_Diff(CMat4D* pMat, const CMat4D* pIn)=0;
/**
\brief subtrai pRight de pLeft e armazena o resultado em pMat. (pMat = pLeft - pRight)
\param pMat matriz que armazenará o resultado da subtração.
\param pLeft matriz que será utilizada na subtração (pMat = pLeft - pRight)
\param pRight matriz que será utilizada na subtração (pMat = pLeft - pRight)
\brief Subtracts pRight from pLeft and stores the result in pMat.
**/
virtual void  VPCALL mat4D_DiffOf(CMat4D* pMat, const CMat4D* pLeft, const CMat4D* pRight)=0;
/**
\brief multiplica pMat por scalare armazena o resultado em pMat. (pMat = pMat * scalar)
\param pMat matriz que será multiplicada e armazenará o resultado
\param scalar númerom que será multiplicador de pMat.
\note Scales the components of pMat by scalar and stores the result in pMat.
**/
virtual void  VPCALL mat4D_Scale(CMat4D* pMat, float scalar)=0;
/**
\brief multiplica pMat por scalare armazena o resultado em pMat. (pMat = pMat * scalar)
\param pOut matriz que  armazenará o resultado
\param pIn matriz que será multiplicada por scalar
\param scalar número que será multiplicador de pMat.
\note Scales the components of pIn by scalar and stores the result in pMat.
**/
virtual void  VPCALL mat4D_ScaleOf(CMat4D* pOut, const CMat4D* pIn, float scalar)=0;
/**
\brief Multiplica pLeft por pRight e armazena o resultado em pLeft. (pLeft = pRight * pLeft)
\param pLeft matriz
\param pRight matriz
\note Multiplies pLeft by pRight and stores the result in pLeft.
**/
virtual void  VPCALL mat4D_Multiply(CMat4D* pLeft, const CMat4D* pRight)=0;
/**
\brief Multiplica pLeft por pRighte armazena o resultado em pOut. (pOut = pRight * pLeft)
\param pLeft matriz
\param pRight matriz
\param pOut Matriz que armazena o resultado da multiplicação
\brief Multiplies pLeft by pRight and stores the result in pLeft.
**/
virtual void  VPCALL mat4D_MultiplyOf(CMat4D* pOut, const CMat4D* pLeft, const CMat4D* pRight)=0;
/**
\brief calcula a matriz transposta de pIn e armazena o resultado em pIn. (pIn = pInt)
\param pIn matriz que se vai calcular a transposta e armazenará o resultado
\note Transposes the matrix pIn stores the result in pIn.
**/
virtual void  VPCALL mat4D_Transpose(CMat4D* pIn)=0;
/**
\brief calcula a matriz transposta de pIn e armazena o resultado em pOut. (pOut = pInt)
\param pIn matriz que se vai calcular a transposta
\param pOut matriz que armazenará o resultado
\note Transposes the matrix pIn stores the result in pOut.
**/
virtual void  VPCALL mat4D_TransposeOf(CMat4D* pOut, const CMat4D* pIn)=0;
/**
\brief Multiplica o vetor pela matriz (pOut = pOut * pIn)
\param pOut Vetor multiplicador
\param pIn Matriz multiplicadora
\note Transforms the 3D vector (as 4D with w = 1.0) pVec by the matrix pMat and stores the result in pVec.
The w component is discarded. If you need to save the w component (for example, a software renderer uses it for perspective correction)
then use 4D vectors with mat4D_VectorMultiply(CVec4D* pOut4D, const CMat4D* pIn).

**/
virtual void  VPCALL mat4D_VectorMultiply(CVec3D* pOut, const CMat4D* pIn)=0;
/**
\brief Multiplica o vetor pela matriz (pOut = pOut * pIn)
\param pOut Vetor que armazenará o resultado
\param pIn Vetor multiplicador
\param pMat Matriz multiplicadora
\note Transforms the 3D vector (as 4D with w = 1.0) pIn by the matrix pMat and stores the result in pOut.
The w component is discarded. If you need to save the w component (for example, a software renderer uses it for perspective correction)
then use 4D vectors with mat4D_VectorMultiply(CVec4D* pOut4D, const CVec4D* pIn4D, const CMat4D* pMat).

**/
virtual void  VPCALL mat4D_VectorMultiplyOf(CVec3D* pOut, const CVec3D* pIn, const CMat4D* pMat)=0;
/**
\brief Multiplica o vetor pela matriz (pOut4D = pOut4D * pIn)
\param pOut4D Vetor multiplicador
\param pIn Matriz multiplicadora
\note Transforms the 4D vector pVec4D by the matrix pMat and stores the result in pVec4D.
**/
virtual void  VPCALL mat4D_VectorMultiply(CVec4D* pOut4D, const CMat4D* pIn)=0;
/**
\brief Multiplica o vetor pela matriz (pOut4D = pIn4D * pMat)
\param pOut4D Vetor que armazenará o resultado
\param pIn4D Vetor multiplicador
\param pMat Matriz multiplicadora
\note Transforms the 4D vector pIn4D by the matrix pMat and stores the result in pOut4D.
**/
virtual void  VPCALL mat4D_VectorMultiplyOf(CVec4D* pOut4D, const CVec4D* pIn4D, const CMat4D* pMat)=0;
/**
\brief rotaciona pMat
\param yaw  parâmetro de rotação no sentido yaw
\param pitch parâmetro de rotação no sentido pitch
\param roll parâmetro de rotação no sentido row
**/
virtual void  VPCALL mat4D_ToRotate(CMat4D* pMat, float yaw, float pitch, float roll)=0;
/**
\brief rotaciona pMat
\param pYawPitchRoll  Vetor que contém os parâmetros de rotação no sentido yaw, pitch, roll de rotação

**/
virtual void  VPCALL mat4D_ToRotateOf(CMat4D* pMat, const CVec3D* pYawPitchRoll)=0;

/**
\brief soma pMat e pIn e guarda o resultado em pMat
\note pMat,pIn devem estar alinhados em 16bytes
\param pMat Matriz que será somada a pIn e guardará o resultado da soma
\param pIn Matria que será somada a pMat
\note Sums pMat and pIn and stores the result in pMat.
\note pMat, pIn must be 16 bytes aligned
**/
virtual void  VPCALL mat4D_AlignedSum(CMat4D* pMat, const CMat4D* pIn)=0;
/**
\brief subtrai pIn de pMat e armazena o resultado em pMat. (pMat = pMat - pIn)
\brief Sums pIn1 and pIn2 and stores the result in pMat.
\note  pMat,pIn1,pIn2 devem estar alinhados em 16bytes.
\param pMat matriz que será utilizada para subtrair pIn e armazenará o resultado da subtração
\param pIn matriz que será subtrída de pMat
\note  pMat,pIn1,pIn2 must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedSumOf(CMat4D* pMat, const CMat4D* pIn1, const CMat4D* pIn2)=0;
/**
\brief subtrai pIn de pMat e armazena o resultado em pMat. (pMat = pMat - pIn)
\brief Subtracts pIn from pMat and stores the result in pMat.
\param pMat matriz que será utilizada para subtrair pIn e armazenará o resultado da subtração
\param pIn matriz que será subtrída de pMat
\note  pMat,pIn devem estar alinhados em 16bytes.
\note pMat,pIn must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedDiff(CMat4D* pMat, const CMat4D* pIn)=0;
/**
\brief subtrai pRight de pLeft e armazena o resultado em pMat. (pMat = pLeft - pRight)
\param pMat matriz que armazenará o resultado da subtração.
\param pLeft matriz que será utilizada na subtração (pMat = pLeft - pRight)
\param pRight matriz que será utilizada na subtração (pMat = pLeft - pRight)
\note Subtracts pRight from pLeft and stores the result in pMat.
\note pMat,pLeft,pRight devem estar alinhados em 16bytes.
\note pMat,pLeft,pRight must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedDiffOf(CMat4D* pMat, const CMat4D* pLeft, const CMat4D* pRight)=0;
/**
\brief multiplica pMat por scalare armazena o resultado em pMat. (pMat = pMat * scalar)
\param pMat matriz que será multiplicada e armazenará o resultado
\param scalar númerom que será multiplicador de pMat.
\note Scales the components of pMat by scalar and stores the result in pMat.
\note pMat deve estar alinhada em 16bytes.
\note pMat must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedScale(CMat4D* pMat, float scalar)=0;
/**
\brief multiplica pMat por scalare armazena o resultado em pMat. (pMat = pMat * scalar)
\param pOut matriz que  armazenará o resultado
\param pIn matriz que será multiplicada por scalar
\param scalar número que será multiplicador de pMat.
\note Scales the components of pIn by scalar and stores the result in pMat.
\note pOut, pIn devem estar alinhados em 16bytes.
\note pOut, pIn must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedScaleOf(CMat4D* pOut, const CMat4D* pIn, float scalar)=0;
/**
\brief Multiplica pLeft por pRight e armazena o resultado em pLeft. (pLeft = pRight * pLeft)
\param pLeft matriz
\param pRight matriz
\note Multiplies pLeft by pRight and stores the result in pLeft.
\note pLeft,pRight devem estar alinhados em 16bytes.
\note pLeft,pRight must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedMultiply(CMat4D* pLeft, const CMat4D* pRight)=0;
/**
\brief Multiplica pLeft por pRighte armazena o resultado em pOut. (pOut = pRight * pLeft)
\param pLeft matriz
\param pRight matriz
\param pOut Matriz que armazena o resultado da multiplicação
\note Multiplies pLeft by pRight and stores the result in pLeft.
\note pOut,pLeft,pRight devem estar alinhados em 16bytes.
\note pOut,pLeft,pRight must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedMultiplyOf(CMat4D* pOut, const CMat4D* pLeft, const CMat4D* pRight)=0;
/**
\brief calcula a matriz transposta de pIn e armazena o resultado em pIn. (pIn = pInt)
\param pIn matriz que se vai calcular a transposta e armazenará o resultado
\brief Transposes the matrix pIn stores the result in pIn.
\note pIn deve estar alinhada em 16bytes.
\note pIn must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedTranspose(CMat4D* pIn)=0;
/**
\brief calcula a matriz transposta de pIn e armazena o resultado em pOut. (pOut = pInt)
\param pIn matriz que se vai calcular a transposta
\param pOut matriz que armazenará o resultado
\note Transposes the matrix pIn stores the result in pOut.
\note pOut, pIn devem estar alinhados em 16bytes.
\note pOut, pIn must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedTransposeOf(CMat4D* pOut, const CMat4D* pIn)=0;
/**
\brief Multiplica o vetor pela matriz (pOut = pOut4D * pIn)
\brief Transforms the 4D vector pVec4D by the matrix pMat and stores the result in pVec4D.
\param pOut Vetor multiplicador e que armazenará o resultado
\param pIn Matriz multiplicadora
\brief Transforms the 3D vector (as 4D with w = 1.0) pVec by the matrix pMat and stores the result in pVec. The w component is discarded. If you need to save the w component (for example, a software renderer uses it for perspective correction) then use 4D vectors with SIMDx86Matrix_Vector4Multiply().
\note pOut, pIn devem estar alinhados em 16bytes.
\note pOut, pIn must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedVectorMultiply(CVec3D* pOut, const CMat4D* pIn)=0;
/**
\brief Multiplica o vetor pela matriz (pOut4D = pIn4D * pMat)
\param pOut Vetor que armazenará o resultado
\param pIn Vetor multiplicador
\param pMat Matriz multiplicadora
\brief Transforms the 3D vector (as 4D with w = 1.0) pIn by the matrix pMat and stores the result in pOut. The w component is discarded. If you need to save the w component (for example, a software renderer uses it for perspective correction) then use 4D vectors with SIMDx86Matrix_Vector4MultiplyOf().
\note pOut, pIn, pMat devem estar alinhados em 16bytes.
\note pOut, pIn, pMat must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedVectorMultiplyOf(CVec3D* pOut, const CVec3D* pIn, const CMat4D* pMat)=0;
/**
\brief Multiplica o vetor pela matriz (pOut4D = pOut4D * pMat)
\brief Transforms the 4D vector pVec4D by the matrix pMat and stores the result in pVec4D.
\param pOut4D Vetor multiplicador e que armazenará o resultado
\param pMat Matriz multiplicadora
\note pOut4D, pMat devem estar alinhados em 16bytes.
\note pOut4D, pMat must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedVectorMultiply(CVec4D* pOut4D, const CMat4D* pMat)=0;
/**
\brief Multiplica o vetor pela matriz (pOut4D = pIn4D * pMat)
\param pOut4D Vetor que armazenará o resultado
\param pIn4D Vetor multiplicador
\param pMat Matriz multiplicadora
\note Transforms the 4D vector pIn4D by the matrix pMat and stores the result in pOut4D.
\note pOut4D, pIn4D, pMat must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedVectorMultiplyOf(CVec4D* pOut4D, const CVec4D* pIn4D, const CMat4D* pMat)=0;
/**
\brief rotaciona pMat
\param yaw  parâmetro de rotação no sentido yaw
\param pitch parâmetro de rotação no sentido pitch
\param roll parâmetro de rotação no sentido row
\note pMat deve estar alinhada em 16bytes.
\note pMat must be 16bytes aligned;
**/
virtual void  VPCALL mat4D_AlignedToRotate(CMat4D* pMat, float yaw, float pitch, float roll)=0;
/**
\brief rotaciona pMat
\param pYawPitchRoll  Vetor que contém os parâmetros de rotação no sentido yaw, pitch, roll de rotação
\note pMat deve estar alinhada em 16bytes.
\note pMat must be 16bytes aligned;
**/

virtual void  VPCALL mat4D_AlignedToRotateOf(CMat4D* pMat, const CVec3D* pYawPitchRoll)=0;

/**
\brief from quaternion to matrix 4D
\param vector [1,0,0,0]
\param q quaternion to convert to matrix
\param [out] mat
**/
virtual void  VPCALL quat_to_mat4x4(sf_m128 q, sf_m128 t, sf_m128 *mat)=0;

/// Compute the product M*v, where M is a 4x4 matrix denoted by an array of 4 sf_m128's, and v is a 4x1 vector.
// If we have SSE 4.1, we can use the dpps (dot product) instruction, _mm_dp_ps intrinsic.
// If we have SSE3, we can repeatedly use haddps to accumulate the result.
virtual sf_m128  VPCALL mat4x4_mul_sse(const sf_m128 *matrix, sf_m128 vector)=0;


virtual sf_m128  VPCALL colmajor_mat4x4_mul_sse(const sf_m128 *matrix, sf_m128 vector)=0;

/// Compute the product M*v, where M is a 3x4 matrix denoted by an array of 4 sf_m128's, and v is a 4x1 vector.
virtual sf_m128  VPCALL mat3x4_mul_sse(const sf_m128 *matrix, sf_m128 vector)=0;

virtual CVec3D  VPCALL mat3x4_mul_vec(const sf_m128 *matrix, sf_m128 vector)=0;
/**
\see http://tptp.cc/mirrors/siyobik.info/instruction/DPPS.html
**/
virtual void  VPCALL mat4x4_mul_dpps(sf_m128 *out, const sf_m128 *m1, const sf_m128 *m2)=0;
virtual void  VPCALL mat4x4_mul_dpps_2(sf_m128 *out, const sf_m128 *m1, const sf_m128 *m2)=0;
virtual void  VPCALL mat4x4_mul_dpps_3(sf_m128 *out, const sf_m128 *m1, const sf_m128 *m2)=0;

virtual void  VPCALL mat4x4_mul_sse(sf_m128 *out, const sf_m128 *m1, const sf_m128 *m2)=0;

virtual void  VPCALL mat3x4_mul_sse(sf_m128 *out, const sf_m128 *m1, const sf_m128 *m2)=0;
/**
\brief Calcure the inverse matrix
\param mat matriz to inverse
\param out inverted matrix
\return return determinant
**/
virtual float  VPCALL mat4x4_inverse(const CMat4D *mat, CMat4D *out)=0;

/// Inverts a 3x4 affine transformation matrix (in row-major format) that only consists of rotation (+possibly mirroring) and translation.
virtual void  VPCALL mat3x4_inverse_orthonormal(sf_m128 *mat, sf_m128 *out)=0;

/**
\brief 	Do one iteration of Newton-Rhapson:
		 e_n = 2*e - x*e^2
\see http://en.wikipedia.org/wiki/Newton%27s_method
\see http://www.sosmath.com/calculus/diff/der07/der07.html
**/
virtual sf_m128  VPCALL  newtonRhapsonRecipStep(sf_m128 recip, sf_m128 estimate)=0;
/**
\see http://en.wikipedia.org/wiki/Newton%27s_method
\see http://www.sosmath.com/calculus/diff/der07/der07.html
**/
virtual sf_m128  VPCALL  newtonRhapsonRecip(sf_m128 recip)=0;

/// Computes the determinant of a 4x4 matrix.
virtual float  VPCALL mat4x4_determinant(const CMat4D *row)=0;

/// Computes the determinant of a 3x4 matrix stored in row-major format. (Treated as a square matrix with last row [0,0,0,1])
virtual float  VPCALL mat3x4_determinant(const sf_m128 *row)=0;

virtual void  VPCALL mat3x4_transpose(const sf_m128 *src, sf_m128 *dst)=0;



//=========================Quaternion========================================
/**
\brief Normaiza o Quaternion e armazena o resultado em pQuat
\param pQuat Quaternion a ser normalizado
\brief Normalizes pQuat and stores it in pQuat.
**/
virtual void VPCALL  quaternion_Normalize(CQuaternion* pQuat)=0;
/**
\brief Normaiza o Quaternion e armazena o resultado em pOut
\param pQuat Quaternion a ser normalizado
\param pOut Quaternion que recebera o resultado
\brief Normalizes pQuat and stores it in pOut.
**/
virtual void VPCALL  quaternion_NormalizeOf(CQuaternion* pOut, const CQuaternion* pQuat)=0;
/**
\brief Multiplica os Quaternions e armazena o resultado em pLeft (pLeft = pLeft * pRight)
\param pLeft Quaternion que será multiplicado e receberá o resultado
\param pRight Quaternion que será multiplicado
\brief Multiplies pLeft by pRight and stores the result in pLeft. The result is not explicitly normalized,
       but will be normal if two normal quaternions are used as inputs.
**/
virtual void VPCALL  quaternion_Multiply(CQuaternion* pLeft, const CQuaternion* pRight)=0;
/**
\if pt_br
\brief Multiplica os Quaternions e armazena o resultado em pLeft (pLeft = pLeft * pRight)
\elseif us_en
\brief Multiplies pLeft by pRight and stores the result in pOut. The result is not explicitly normalized,
       but will be normal if two normal quaternions are used as inputs.
\endif
\param pLeft Quaternion que será multiplicado
\param pRight Quaternion que será multiplicado
\param pOut Quaternion que receberá o resultado
**/
virtual void VPCALL  quaternion_MultiplyOf(CQuaternion* pOut, const CQuaternion* pLeft, const CQuaternion* pRight)=0;

#ifdef ANDROID
inline void quat_mul_quat_asm(const void *q1, const void *q2, void *out)
{
/*	return Quat(x*r.w + y*r.z - z*r.y + w*r.x,
	           -x*r.z + y*r.w + z*r.x + w*r.y,
	            x*r.y - y*r.x + z*r.w + w*r.z,
	           -x*r.x - y*r.y - z*r.z + w*r.w); */
#ifdef _DEBUG
	SMF_ASSERT(IS16ALIGNED(out));
	SMF_ASSERT(IS16ALIGNED(q1));
	SMF_ASSERT(IS16ALIGNED(q2));
	SMF_ASSERT(IS16ALIGNED(sx));
	SMF_ASSERT(IS16ALIGNED(sy));
	SMF_ASSERT(IS16ALIGNED(sz));
#endif
	///\todo 128-bit aligned loads: [%1,:128]
	asm(
		"\t vld1.32 {d0, d1}, [%1]\n"    // q0 = quat1.xyzw
		"\t vmov.i32 d12, #0\n"          // q6.lo = 0
		"\t vmov.i32 d13, #0x80000000\n" // q6.hi = [- - + +] = 'signy'
		"\t vld1.32 {d8, d9}, [%2]\n"    // q4 = quat2.xyzw [%2]
		"\t vdup.32 q1, d0[1]\n"         // q1 = q0[1] = quat1.yyyy = 'Y'
		"\t vdup.32 q2, d1[0]\n"         // q2 = q0[2] = quat1.zzzz = 'Z'
		"\t vshl.i64 d10, d13, #32\n"    // q5.lo = q6.hi = [- +]
		"\t vdup.32 q3, d1[1]\n"         // q3 = q0[3] = quat1.wwww = 'W'
		"\t vdup.32 q0, d0[0]\n"         // q0 = q0[0] = quat1.xxxx = 'X'
		"\t vmov d11, d10\n"             // q5.hi = q5.lo = [- + - +] = 'signx'
		"\t vmov d15, d10\n"             // q7.hi = q5.lo = [- +]
		"\t vshr.u64 d14, d10, #32\n"    // q7.lo = q5.lo = [- + + -] = 'signz'

		"\t vmov d18, d9\n"              // q9.lo = q4.hi
		"\t vmov d19, d8\n"              // q9.hi = q4.lo, q9 = quat2.zwxy = 'q2 for Y'

		"\t veor q0, q0, q5\n"           // q0 = X*signx = [-x x -x x]
		"\t veor q1, q1, q6\n"           // q1 = Y*signy = [-y -y y y]
		"\t veor q2, q2, q7\n"           // q2 = Z*signz = [-z z z -z]

		"\t vrev64.32 q10, q9\n"         // q10 = quat2.wzyx = 'q2 for X'

		"\t vmul.f32 q0, q0, q10\n"      // q0 = X*signx * quat2
		"\t vmul.f32 q11, q1, q9\n"       // q0 += Y*signy * quat2
		"\t vrev64.32 q8, q4\n"          // q8 = quat2.yxwz  = 'q2 for Z'
		"\t vmla.f32 q0, q2, q8\n"       // q0 += Z*signz * quat2
		"\t vmla.f32 q11, q3, q4\n"       // q0 += W       * quat2

		"\t vadd.f32 q0, q0, q11\n"
		"\t vst1.32	{d0, d1}, [%0]\n"    // store output
	: /* no outputs by value */
	: [out]"r"(out), [quat1]"r"(q1), [quat2]"r"(q2)
	: "memory", "q11", "q10", "q9", "q8", "q7", "q6", "q5", "q4", "q3", "q2", "q1", "q0");
}
#endif

/// Converts a quaternion to a row-major matrix.
/// From http://renderfeather.googlecode.com/hg-history/034a1900d6e8b6c92440382658d2b01fc732c5de/Doc/optimized%20Matrix%20quaternion%20conversion.pdf
virtual void VPCALL quat_to_mat3x4(sf_m128 q, sf_m128 t, sf_m128 *m)=0;
virtual sf_m128 VPCALL quat_transform_vec4(sf_m128 quat, sf_m128 vec)=0;
virtual sf_m128 VPCALL quat_mul_quat(sf_m128 q1, sf_m128 q2)=0;
virtual sf_m128 VPCALL quat_div_quat(sf_m128 q1, sf_m128 q2)=0;

//=============PLANE=================================

/**
\brief constrói um plano através de três pontos. O Plano não é normalizado
\param pOut O Plano que será criado
\param pA primeiro ponto
\param pB segundo ponto
\param pC terceiro ponto
\brief This function constructs a plane from three points. The plane is not normalized.
**/
virtual void   VPCALL  plane_FromPoints(CPlane* pOut, const CVec3D* pA, const CVec3D* pB,const CVec3D* pC)=0;
/**
\brief retorna a distância sem sinal(módulo), entre um ponto e um plano
\param pPlane Plano a seu utilizado no cálculo
\param pPoint ponto a ser utilizado no cálculo
\brief This function returns the unsigned distance between a point and a plane.
**/
virtual float  VPCALL  plane_DistToPoint(const CPlane* pPlane, const CVec3D* pPoint)=0;
/**
\brief calcula o produto escalar entre o plano e um ponto com w=1.0f. Util para classificar
       um ponto em relação ao plano: um valor menor que zero significa que o ponto está atraz do plano
	   e um valor maoir que zero significa que o ponto está na frente do plano
\param pPlane Plano que será utilizado no cálculo
\param pVec Ponto que será utilizado no cálculo
\brief This function returns the dot product between a plane and a point with w = 1.0.
This is useful for classifying a point in relation to a plane: a value of less than zero
implies the point it behind the plane, a value of zero implies the point is on the plane,
and value of greater than zero implies the point is in front of the plane.

**/
virtual float  VPCALL  plane_Dot(const CPlane* pPlane, const CVec3D* pVec)=0;
/**
\brief retorna o produto escalar entre o plano e um ponto 4D
\param pPlane plano que será utilizado no cálculo do produto
\param pVec4 ponto de 4 dimensões que será utilizado no calculo do produto
\brief This function returns the dot product between a plane and a 4D point.
This is useful for classifying a point in relation to a plane: a value of less than zero implies the point it behind the plane, a value of zero implies the point is on the plane, and value of greater than zero implies the point is in front of the plane.

**/
virtual float  VPCALL  plane_Dot4(const CPlane* pPlane, const CVec4D* pVec4)=0;
/**
\brief retorna o produoto escalar entre o plano e outra normal. representa o cosseno do angulo entre os dois, se ambos estiverem normalizados
\param pPlane Plano que será utilizado para calcular o produto
\param pVec normal do segundo plano
\return O resultado do produto
\brief This functions returns the dot product between the plane's normal and another normal.
This represents the cosine of the angle between the two, if both are normalized.

**/
virtual float  VPCALL  plane_DotNormal(const CPlane* pPlane, const CVec3D* pVec)=0;
/**
\brief calcula o produto escalar entre dois planos normalizados.
\param pPlane1 Plano que será multiplicado
\param pPlane2 Plano que será multiplicado
\return o resultado do produto
\brief This functions returns the dot product between the two planes' normals.
This represents the cosine of the dihedral angle between the two, if both are normalized.
**/
virtual float  VPCALL  plane_DotPlane(const CPlane* pPlane1, const CPlane* pPlane2)=0;
/**
\brief Normaliza o plano pOut e armazena o resultado nele mesmo
\param pOut plano que será normalizado e armazenará o resultado
\brief This function normalizes the plane such that the magnitude of its normal is 1 and modifies 'd' component appropriately.

**/
virtual void   VPCALL  plane_Normalize(CPlane* pOut)=0;
/**
\brief Normaliza o plano e armazena o resultado em pOut
\param pOut plano que armazenará o resultado
\param pIn Plano que será normalizado
\brief This function stores the normalized plane pIn into pOut, while pIn is unaffected.
**/
virtual void  VPCALL  plane_NormalizeOf(CPlane* pOut, CPlane* pIn)=0;

virtual	bool  VPCALL  intersectLineAABB(const CAABBox &box, const CVec4D &rayPos, const CVec4D &rayDir, float tNear, float tFar)=0;


public:

	cpuid_t	cpuid;

};

// pointer to SIMD processor
extern CSIMDProcessor *SIMDProcessor;
} //end MATH
} //end SMF
#endif /* !__SMF_MATH_SIMD_H__ */
