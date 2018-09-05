
#include <io.h>
#include <string>
#include <vector>
#include <queue>
#include <list>
#include <iostream>
#include "SMF.h"
#include "math/all.h"

using namespace std;
using std::set_terminate;
using std::set_unexpected;
using std::set_new_handler;

#undef printf


#define SAMPLES 100


int i = 0;
long base = 0;
long tick = 0;
long ticks[SAMPLES];

int Duration(int sz = SAMPLES)
{
long nclocks = 0;
	for (int i = 0 ; i < sz; i++) {
	if (!nclocks || ticks[i] < nclocks)
	nclocks = ticks[i];
	}

return int(nclocks - base);

}

void report(char* format, ...)
{
va_list marker;
char buf[500];
vsprintf(buf, format, va_start(marker, format));
//puts(buf);
SMF::Debug::debug(SMF::Debug::math,__FUNCTION__) << buf << endl;
// OutputDebugString(buf); OutputDebugString("\n");

}

#define START_MEASUREMENTS \
SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_TIME_CRITICAL); \
for (i = 0 ; i < SAMPLES; i++) {

#define END_MEASUREMENTS \
ticks[i] = SMF::MATH::end_ClockCount - SMF::MATH::start_ClockCount; \
} \
SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_NORMAL); \
report("Duration for %s:\t%i", testname, Duration());

// Offset for mat[i][j], w is an row width, t == 1 for transposed access e t==0 for normal access.
#define mi(w, t, i, j) 4 * ((i * w + j) * (1-t) + (j * w + i) * t)

// Load & multiply.
#define flm(k, i, j, m, n, a, b) \
__asm fld dword ptr [edx + mi(m, a, i, k)] \
__asm fmul dword ptr [ecx + mi(n, b, k, j)]

// Load, multiply & add.
#define flma(k, i, j, m, n, a, b) flm(k, i, j, m, n, a, b) __asm faddp ST(1), ST(0)

#define e3(i, j, l, m, n, a, b) \
flm (0, i, j, m, n, a, b) \
flma(1, i, j, m, n, a, b) \
flma(2, i, j, m, n, a, b) \
__asm fstp dword ptr [eax + mi(l, 0, i, j)]

void PII_Mult_3x3_3x3(float *src1, float *src2, float *dst)
{
StartRecordTime;
__asm mov edx, DWORD PTR src1
__asm mov ecx, DWORD PTR src2
__asm mov eax, DWORD PTR dst
e3(0, 0, 3, 3, 3, 0, 0) e3(0, 1, 3, 3, 3, 0, 0) e3(0, 2, 3, 3, 3, 0, 0)
e3(1, 0, 3, 3, 3, 0, 0) e3(1, 1, 3, 3, 3, 0, 0) e3(1, 2, 3, 3, 3, 0, 0)
e3(2, 0, 3, 3, 3, 0, 0) e3(2, 1, 3, 3, 3, 0, 0) e3(2, 2, 3, 3, 3, 0, 0)
StopRecordTime;
}

#define e4(i, j, l, m, n, a, b) \
flm(0, i, j, m, n, a, b) \
flm(1, i, j, m, n, a, b) \
flm(2, i, j, m, n, a, b) \
flm(3, i, j, m, n, a, b) \
__asm faddp st(1), st(0) \
__asm fxch st(2) \
__asm faddp st(1), st(0) \
__asm faddp st(1), st(0) \
__asm fstp dword ptr [eax + mi(l, 0, i, j)]

void PII_Mult00_4x4_4x4(float *src1, float *src2, float *dst)
{
StartRecordTime;
__asm mov edx, DWORD PTR src1
__asm mov ecx, DWORD PTR src2
__asm mov eax, DWORD PTR dst
e4(0, 0, 4, 4, 4, 0, 0)
e4(0, 1, 4, 4, 4, 0, 0)
e4(0, 2, 4, 4, 4, 0, 0)
e4(0, 3, 4, 4, 4, 0, 0)
e4(1, 0, 4, 4, 4, 0, 0)
e4(1, 1, 4, 4, 4, 0, 0)
e4(1, 2, 4, 4, 4, 0, 0)
e4(1, 3, 4, 4, 4, 0, 0)
e4(2, 0, 4, 4, 4, 0, 0)
e4(2, 1, 4, 4, 4, 0, 0)
e4(2, 2, 4, 4, 4, 0, 0)
e4(2, 3, 4, 4, 4, 0, 0)
e4(3, 0, 4, 4, 4, 0, 0)
e4(3, 1, 4, 4, 4, 0, 0)
e4(3, 2, 4, 4, 4, 0, 0)
e4(3, 3, 4, 4, 4, 0, 0)
StopRecordTime;
}

void PII_Mult00_4x4_4x1(float *src1, float *src2, float *dst)

{
StartRecordTime;
__asm mov edx, DWORD PTR src1
__asm mov ecx, DWORD PTR src2
__asm mov eax, DWORD PTR dst
e4(0, 0, 1, 4, 1, 0, 0)
e4(1, 0, 1, 4, 1, 0, 0)
e4(2, 0, 1, 4, 1, 0, 0)
e4(3, 0, 1, 4, 1, 0, 0)
StopRecordTime;
}

#define e6(i, j, l, m, n, a, b) \
flm(0, i, j, m, n, a, b) \
flm(1, i, j, m, n, a, b) \
flm(2, i, j, m, n, a, b) \
flm(3, i, j, m, n, a, b) \
flm(4, i, j, m, n, a, b) \
flm(5, i, j, m, n, a, b) \
__asm faddp st(1), st(0) \
__asm fxch st(2) \
__asm faddp st(1), st(0) \
__asm faddp st(1), st(0) \
__asm fxch st(2) \
__asm faddp st(1), st(0) \
__asm faddp st(1), st(0) \
__asm fstp dword ptr [eax + mi(l, 0, i, j)]

void PII_Mult00_6x6_6x6(float *src1, float *src2, float *dst)
{
StartRecordTime;
__asm mov edx, DWORD PTR src1
__asm mov ecx, DWORD PTR src2
__asm mov eax, DWORD PTR dst
e6(0, 0, 6, 6, 6, 0, 0)
e6(0, 1, 6, 6, 6, 0, 0)
e6(0, 2, 6, 6, 6, 0, 0)
e6(0, 3, 6, 6, 6, 0, 0)
e6(0, 4, 6, 6, 6, 0, 0)
e6(0, 5, 6, 6, 6, 0, 0)
e6(1, 0, 6, 6, 6, 0, 0)
e6(1, 1, 6, 6, 6, 0, 0)
e6(1, 2, 6, 6, 6, 0, 0)
e6(1, 3, 6, 6, 6, 0, 0)
e6(1, 4, 6, 6, 6, 0, 0)
e6(1, 5, 6, 6, 6, 0, 0)
e6(2, 0, 6, 6, 6, 0, 0)
e6(2, 1, 6, 6, 6, 0, 0)
e6(2, 2, 6, 6, 6, 0, 0)
e6(2, 3, 6, 6, 6, 0, 0)
e6(2, 4, 6, 6, 6, 0, 0)
e6(2, 5, 6, 6, 6, 0, 0)
e6(3, 0, 6, 6, 6, 0, 0)
e6(3, 1, 6, 6, 6, 0, 0)
e6(3, 2, 6, 6, 6, 0, 0)
e6(3, 3, 6, 6, 6, 0, 0)
e6(3, 4, 6, 6, 6, 0, 0)
e6(3, 5, 6, 6, 6, 0, 0)
e6(4, 0, 6, 6, 6, 0, 0)
e6(4, 1, 6, 6, 6, 0, 0)
e6(4, 2, 6, 6, 6, 0, 0)
e6(4, 3, 6, 6, 6, 0, 0)
e6(4, 4, 6, 6, 6, 0, 0)
e6(4, 5, 6, 6, 6, 0, 0)
e6(5, 0, 6, 6, 6, 0, 0)
e6(5, 1, 6, 6, 6, 0, 0)
e6(5, 2, 6, 6, 6, 0, 0)
e6(5, 3, 6, 6, 6, 0, 0)
e6(5, 4, 6, 6, 6, 0, 0)
e6(5, 5, 6, 6, 6, 0, 0)
StopRecordTime;
}

void PII_Mult00_6x6_6x1(float *src1, float *src2, float *dst)
{
StartRecordTime;
__asm mov edx, DWORD PTR src1
__asm mov ecx, DWORD PTR src2
__asm mov eax, DWORD PTR dst
e6(0, 0, 1, 6, 1, 0, 0)
e6(1, 0, 1, 6, 1, 0, 0)
e6(2, 0, 1, 6, 1, 0, 0)
e6(3, 0, 1, 6, 1, 0, 0)
e6(4, 0, 1, 6, 1, 0, 0)
e6(5, 0, 1, 6, 1, 0, 0)
StopRecordTime;
}

void PII_Mult00_3x3_3x1(float *src1, float *src2, float *dst)
{
StartRecordTime;
__asm {
mov edx, dword ptr src1
mov ecx, dword ptr src2
mov eax, dword ptr dst
fld dword ptr [ecx]
fmul dword ptr [edx+24]
fld dword ptr [ecx]
fmul dword ptr [edx+12]
fld dword ptr [ecx]
fmul dword ptr [edx]
fld dword ptr [ecx+4]
fmul dword ptr [edx+4]
fld dword ptr [ecx+4]
fmul dword ptr [edx+16]
fld dword ptr [ecx+4]
fmul dword ptr [edx+28]
fxch ST(2)
faddp ST(3),ST
faddp ST(3),ST
faddp ST(3),ST
fld dword ptr [ecx+8]
fmul dword ptr [edx+8]
fld dword ptr [ecx+8]
fmul dword ptr [edx+20]
fld dword ptr [ecx+8]
fmul dword ptr [edx+32]
fxch ST(2)
faddp ST(3),ST
faddp ST(3),ST
faddp ST(3),ST
fstp dword ptr [eax]
fstp dword ptr [eax+4]
fstp dword ptr [eax+8]
}
StopRecordTime;
}

__declspec(naked) void PIII_Mult00_3x3_3x3(float *src1, float *src2, float *dst)
{
StartRecordTime;
__asm {
mov ecx, dword ptr [esp+8] ; src2
mov edx, dword ptr [esp+4] ; src1
mov eax, dword ptr [esp+0Ch] ; dst
movss xmm2, dword ptr [ecx+32]
movhps xmm2, qword ptr [ecx+24]
movss xmm3, dword ptr [edx]
movss xmm4, dword ptr [edx+4]
movss xmm0, dword ptr [ecx]
movhps xmm0, qword ptr [ecx+4]
shufps xmm2, xmm2, 0x36
shufps xmm3, xmm3, 0
movss xmm1, dword ptr [ecx+12]
movhps xmm1, qword ptr [ecx+16]
shufps xmm4, xmm4, 0
mulps xmm3, xmm0
movss xmm5, dword ptr [edx+8]
movss xmm6, dword ptr [edx+12]
mulps xmm4, xmm1
shufps xmm5, xmm5, 0
mulps xmm5, xmm2
shufps xmm6, xmm6, 0
mulps xmm6, xmm0
addps xmm3, xmm4
movss xmm7, dword ptr [edx+16]
movss xmm4, dword ptr [edx+28]
shufps xmm7, xmm7, 0
addps xmm3, xmm5
mulps xmm7, xmm1
shufps xmm4, xmm4, 0
movss xmm5, dword ptr [edx+20]
shufps xmm5, xmm5, 0
mulps xmm4, xmm1
mulps xmm5, xmm2
addps xmm6, xmm7
movss xmm1, dword ptr [edx+24]
movss dword ptr [eax] , xmm3
movhps qword ptr [eax+4], xmm3
addps xmm6, xmm5
shufps xmm1, xmm1, 0
movss xmm5, dword ptr [edx+32]
mulps xmm1, xmm0
shufps xmm5, xmm5, 0
movss dword ptr [eax+12], xmm6
mulps xmm5, xmm2
addps xmm1, xmm4
movhps qword ptr [eax+16], xmm6
addps xmm1, xmm5
shufps xmm1, xmm1, 0x8F
movhps qword ptr [eax+24], xmm1
movss dword ptr [eax+32], xmm1
}
StopRecordTime;
__asm ret
}

__declspec(naked) void PIII_Mult00_4x4_4x4(float *src1, float *src2, float *dst)
{
StartRecordTime;
__asm {
mov edx, dword ptr [esp+4] ; src1
mov eax, dword ptr [esp+0Ch] ; dst
mov ecx, dword ptr [esp+8] ; src2
movss xmm0, dword ptr [edx]
movaps xmm1, xmmword ptr [ecx]
shufps xmm0, xmm0, 0
movss xmm2, dword ptr [edx+4]
mulps xmm0, xmm1
shufps xmm2, xmm2, 0
movaps xmm3, xmmword ptr [ecx+10h]
movss xmm7, dword ptr [edx+8]
mulps xmm2, xmm3
shufps xmm7, xmm7, 0
addps xmm0, xmm2
movaps xmm4, xmmword ptr [ecx+20h]
movss xmm2, dword ptr [edx+0Ch]
mulps xmm7, xmm4
shufps xmm2, xmm2, 0
addps xmm0, xmm7
movaps xmm5, xmmword ptr [ecx+30h]
movss xmm6, dword ptr [edx+10h]
mulps xmm2, xmm5
movss xmm7, dword ptr [edx+14h]
shufps xmm6, xmm6, 0
addps xmm0, xmm2
shufps xmm7, xmm7, 0
movlps qword ptr [eax], xmm0
movhps qword ptr [eax+8], xmm0
mulps xmm7, xmm3
movss xmm0, dword ptr [edx+18h]
mulps xmm6, xmm1
shufps xmm0, xmm0, 0
addps xmm6, xmm7
mulps xmm0, xmm4
movss xmm2, dword ptr [edx+24h]
addps xmm6, xmm0
movss xmm0, dword ptr [edx+1Ch]
movss xmm7, dword ptr [edx+20h]
shufps xmm0, xmm0, 0
shufps xmm7, xmm7, 0
mulps xmm0, xmm5
mulps xmm7, xmm1
addps xmm6, xmm0
shufps xmm2, xmm2, 0
movlps qword ptr [eax+10h], xmm6
movhps qword ptr [eax+18h], xmm6
mulps xmm2, xmm3
movss xmm6, dword ptr [edx+28h]
addps xmm7, xmm2
shufps xmm6, xmm6, 0
movss xmm2, dword ptr [edx+2Ch]
mulps xmm6, xmm4
shufps xmm2, xmm2, 0
addps xmm7, xmm6
mulps xmm2, xmm5
movss xmm0, dword ptr [edx+34h]
addps xmm7, xmm2
shufps xmm0, xmm0, 0
movlps qword ptr [eax+20h], xmm7
movss xmm2, dword ptr [edx+30h]
movhps qword ptr [eax+28h], xmm7
mulps xmm0, xmm3
shufps xmm2, xmm2, 0
movss xmm6, dword ptr [edx+38h]
mulps xmm2, xmm1
shufps xmm6, xmm6, 0
addps xmm2, xmm0
mulps xmm6, xmm4
movss xmm7, dword ptr [edx+3Ch]
shufps xmm7, xmm7, 0
addps xmm2, xmm6
mulps xmm7, xmm5
addps xmm2, xmm7
movaps xmmword ptr [eax+30h], xmm2
}
StopRecordTime;
__asm ret
}

__declspec(naked) void PIII_Mult00_6x6_6x6(float *src1, float *src2, float *dst)
{
StartRecordTime;
__asm {
mov ecx, dword ptr [esp+8] ; src2
movlps xmm3, qword ptr [ecx+72]
mov edx, dword ptr [esp+4] ; src1
// Loading first 4 columns (upper 4 rows) of src2.
movaps xmm0, xmmword ptr [ecx]
movlps xmm1, qword ptr [ecx+24]
movhps xmm1, qword ptr [ecx+32]
movaps xmm2, xmmword ptr [ecx+48]
movhps xmm3, qword ptr [ecx+80]
// Calculating first 4 elements in the first row of the destination matrix.
movss xmm4, dword ptr [edx]
movss xmm5, dword ptr [edx+4]
mov eax, dword ptr [esp+0Ch] ; dst
shufps xmm4, xmm4, 0
movss xmm6, dword ptr [edx+8]
shufps xmm5, xmm5, 0
movss xmm7, dword ptr [edx+12]
mulps xmm4, xmm0
shufps xmm6, xmm6, 0
shufps xmm7, xmm7, 0
mulps xmm5, xmm1
mulps xmm6, xmm2
addps xmm5, xmm4
mulps xmm7, xmm3
addps xmm6, xmm5
addps xmm7, xmm6
movaps xmmword ptr [eax], xmm7
// Calculating first 4 elements in the second row of the destination matrix.
movss xmm4, dword ptr [edx+24]
shufps xmm4, xmm4, 0
mulps xmm4, xmm0
movss xmm5, dword ptr [edx+28]
shufps xmm5, xmm5, 0
mulps xmm5, xmm1
movss xmm6, dword ptr [edx+32]
shufps xmm6, xmm6, 0
movss xmm7, dword ptr [edx+36]
shufps xmm7, xmm7, 0
mulps xmm6, xmm2
mulps xmm7, xmm3
addps xmm7, xmm6
addps xmm5, xmm4
addps xmm7, xmm5
// Calculating first 4 elements in the third row of the destination matrix.
movss xmm4, dword ptr [edx+48]
movss xmm5, dword ptr [edx+52]
movlps qword ptr [eax+24], xmm7 ; save 2nd
movhps qword ptr [eax+32], xmm7 ; row
movss xmm6, dword ptr [edx+56]
movss xmm7, dword ptr [edx+60]
shufps xmm4, xmm4, 0
shufps xmm5, xmm5, 0
shufps xmm6, xmm6, 0
shufps xmm7, xmm7, 0
mulps xmm4, xmm0
mulps xmm5, xmm1
mulps xmm6, xmm2
mulps xmm7, xmm3
addps xmm5, xmm4
addps xmm7, xmm6
addps xmm7, xmm5
movaps xmmword ptr [eax+48], xmm7

// Calculating first 4 elements in the fourth row of the destination matrix.

movss xmm4, dword ptr [edx+72]
movss xmm5, dword ptr [edx+76]
movss xmm6, dword ptr [edx+80]
movss xmm7, dword ptr [edx+84]
shufps xmm4, xmm4, 0
shufps xmm5, xmm5, 0
shufps xmm6, xmm6, 0
shufps xmm7, xmm7, 0
mulps xmm4, xmm0
mulps xmm5, xmm1
mulps xmm6, xmm2
mulps xmm7, xmm3
addps xmm4, xmm5
addps xmm6, xmm4
addps xmm7, xmm6
movlps qword ptr [eax+72], xmm7
movhps qword ptr [eax+80], xmm7

// Calculating first 4 elements in the fifth row of the destination matrix.

movss xmm4, dword ptr [edx+96]
movss xmm5, dword ptr [edx+100]
movss xmm6, dword ptr [edx+104]
movss xmm7, dword ptr [edx+108]
shufps xmm4, xmm4, 0
shufps xmm5, xmm5, 0
shufps xmm6, xmm6, 0
shufps xmm7, xmm7, 0
mulps xmm4, xmm0
mulps xmm5, xmm1
mulps xmm6, xmm2
mulps xmm7, xmm3
addps xmm5, xmm4
addps xmm7, xmm6
addps xmm7, xmm5
movaps xmmword ptr [eax+96], xmm7

// Calculating first 4 elements in the sixth row of the destination matrix.

movss xmm4, dword ptr [edx+120]
movss xmm5, dword ptr [edx+124]
movss xmm6, dword ptr [edx+128]
movss xmm7, dword ptr [edx+132]
shufps xmm4, xmm4, 0
shufps xmm5, xmm5, 0
shufps xmm6, xmm6, 0
shufps xmm7, xmm7, 0
mulps xmm4, xmm0
mulps xmm5, xmm1
mulps xmm6, xmm2
mulps xmm7, xmm3
addps xmm4, xmm5
addps xmm6, xmm4
addps xmm7, xmm6
movhps qword ptr [eax+128], xmm7
movlps qword ptr [eax+120], xmm7

// Loading first 4 columns (lower 2 rows) of src2.

movlps xmm0, qword ptr [ecx+96]
movhps xmm0, qword ptr [ecx+104]
movlps xmm1, qword ptr [ecx+120]
movhps xmm1, qword ptr [ecx+128]

// Calculating first 4 elements in the first row of the destination matrix.

movss xmm2, dword ptr [edx+16]
shufps xmm2, xmm2, 0
movss xmm4, dword ptr [edx+40]
movss xmm3, dword ptr [edx+20]
movss xmm5, dword ptr [edx+44]
movaps xmm6, xmmword ptr [eax]
movlps xmm7, qword ptr [eax+24]
shufps xmm3, xmm3, 0
shufps xmm5, xmm5, 0
movhps xmm7, qword ptr [eax+32]
shufps xmm4, xmm4, 0
mulps xmm5, xmm1
mulps xmm2, xmm0
mulps xmm3, xmm1
mulps xmm4, xmm0
addps xmm6, xmm2
addps xmm7, xmm4
addps xmm7, xmm5
addps xmm6, xmm3
movlps qword ptr [eax+24], xmm7
movaps xmmword ptr [eax], xmm6
movhps qword ptr [eax+32], xmm7

// Calculating first 4 elements in the third row of the destination matrix.

movss xmm2, dword ptr [edx+64]
movss xmm4, dword ptr [edx+88]
movss xmm5, dword ptr [edx+92]
movss xmm3, dword ptr [edx+68]
movaps xmm6, xmmword ptr [eax+48]
movlps xmm7, qword ptr [eax+72]
movhps xmm7, qword ptr [eax+80]
shufps xmm2, xmm2, 0
shufps xmm4, xmm4, 0
shufps xmm5, xmm5, 0
shufps xmm3, xmm3, 0
mulps xmm2, xmm0
mulps xmm4, xmm0
mulps xmm5, xmm1
mulps xmm3, xmm1
addps xmm6, xmm2
addps xmm6, xmm3
addps xmm7, xmm4
addps xmm7, xmm5
movlps qword ptr [eax+72], xmm7
movaps xmmword ptr [eax+48], xmm6
movhps qword ptr [eax+80], xmm7

// Calculating first 4 elements in the fifth row of the destination matrix.

movss xmm2, dword ptr [edx+112]
movss xmm3, dword ptr [edx+116]
movaps xmm6, xmmword ptr [eax+96]
shufps xmm2, xmm2, 0
shufps xmm3, xmm3, 0
mulps xmm2, xmm0
mulps xmm3, xmm1
addps xmm6, xmm2
addps xmm6, xmm3
movaps xmmword ptr [eax+96], xmm6

// Calculating first 4 elements in the sixth row of the destination matrix.

movss xmm4, dword ptr [edx+136]
movss xmm5, dword ptr [edx+140]
movhps xmm7, qword ptr [eax+128]
movlps xmm7, qword ptr [eax+120]
shufps xmm4, xmm4, 0
shufps xmm5, xmm5, 0
mulps xmm4, xmm0
mulps xmm5, xmm1
addps xmm7, xmm4
addps xmm7, xmm5

// Calculating last 2 columns of the destination matrix.

movlps xmm0, qword ptr [ecx+16]
movhps xmm0, qword ptr [ecx+40]
movhps qword ptr [eax+128], xmm7
movlps qword ptr [eax+120], xmm7
movlps xmm2, qword ptr [ecx+64]
movhps xmm2, qword ptr [ecx+88]
movaps xmm3, xmm2
shufps xmm3, xmm3, 4Eh
movlps xmm4, qword ptr [ecx+112]
movhps xmm4, qword ptr [ecx+136]
movaps xmm5, xmm4
shufps xmm5, xmm5, 4Eh
movlps xmm6, qword ptr [edx]
movhps xmm6, qword ptr [edx+24]
movaps xmm7, xmm6
shufps xmm7, xmm7, 0F0h
mulps xmm7, xmm0
shufps xmm6, xmm6, 0A5h
movaps xmm1, xmm0
shufps xmm1, xmm1, 4Eh
mulps xmm1, xmm6
addps xmm7, xmm1
movlps xmm6, qword ptr [edx+8]
movhps xmm6, qword ptr [edx+32]
movaps xmm1, xmm6
shufps xmm1, xmm1, 0F0h
shufps xmm6, xmm6, 0A5h
mulps xmm1, xmm2
mulps xmm6, xmm3
addps xmm7, xmm1
addps xmm7, xmm6
movhps xmm6, qword ptr [edx+40]
movlps xmm6, qword ptr [edx+16]
movaps xmm1, xmm6
shufps xmm1, xmm1, 0F0h
shufps xmm6, xmm6, 0A5h
mulps xmm1, xmm4
mulps xmm6, xmm5
addps xmm7, xmm1
addps xmm7, xmm6
movlps qword ptr [eax+16], xmm7
movhps qword ptr [eax+40], xmm7
movlps xmm6, qword ptr [edx+48]
movhps xmm6, qword ptr [edx+72]
movaps xmm7, xmm6
shufps xmm7, xmm7, 0F0h
mulps xmm7, xmm0
shufps xmm6, xmm6, 0A5h
movaps xmm1, xmm0
shufps xmm1, xmm1, 4Eh
mulps xmm1, xmm6
addps xmm7, xmm1
movhps xmm6, qword ptr [edx+80]
movlps xmm6, qword ptr [edx+56]
movaps xmm1, xmm6
shufps xmm1, xmm1, 0F0h
shufps xmm6, xmm6, 0A5h
mulps xmm1, xmm2
mulps xmm6, xmm3
addps xmm7, xmm1
addps xmm7, xmm6
movlps xmm6, qword ptr [edx+64]
movhps xmm6, qword ptr [edx+88]
movaps xmm1, xmm6
shufps xmm1, xmm1, 0F0h
shufps xmm6, xmm6, 0A5h
mulps xmm1, xmm4
mulps xmm6, xmm5
addps xmm7, xmm1
addps xmm7, xmm6
movlps qword ptr [eax+64], xmm7
movhps qword ptr [eax+88], xmm7
movlps xmm6, qword ptr [edx+96]
movhps xmm6, qword ptr [edx+120]
movaps xmm7, xmm6
shufps xmm7, xmm7, 0F0h
mulps xmm7, xmm0
shufps xmm6, xmm6, 0A5h
movaps xmm1, xmm0
shufps xmm1, xmm1, 4Eh
mulps xmm1, xmm6
addps xmm7, xmm1
movlps xmm6, qword ptr [edx+104]
movhps xmm6, qword ptr [edx+128]
movaps xmm1, xmm6
shufps xmm1, xmm1, 0F0h
shufps xmm6, xmm6, 0A5h
mulps xmm1, xmm2
mulps xmm6, xmm3
addps xmm7, xmm1
addps xmm7, xmm6
movlps xmm6, qword ptr [edx+112]
movhps xmm6, qword ptr [edx+136]
movaps xmm1, xmm6
shufps xmm1, xmm1, 0F0h
shufps xmm6, xmm6, 0A5h
mulps xmm1, xmm4
mulps xmm6, xmm5
addps xmm7, xmm1
addps xmm7, xmm6
movlps qword ptr [eax+112], xmm7
movhps qword ptr [eax+136], xmm7
}
StopRecordTime;
__asm ret
}

__declspec(naked) void PIII_Mult00_3x3_3x1(float *src1, float *src2, float *dst)
{
StartRecordTime;
__asm {
mov edx, dword ptr [esp+4] ; src1
mov ecx, dword ptr [esp+8] ; src2
movss xmm1, dword ptr [edx]
mov eax, dword ptr [esp+0Ch] ; dst
movhps xmm1, qword ptr [edx+4]
movaps xmm5, xmm1
movss xmm3, dword ptr [edx+12]
movhps xmm3, qword ptr [edx+24]
movss xmm4, dword ptr [ecx]
shufps xmm5, xmm3, 128
movlps xmm0, qword ptr [edx+16]
shufps xmm4, xmm4, 0
movhps xmm0, qword ptr [edx+28]
shufps xmm1, xmm0, 219
movss xmm2, dword ptr [ecx+4]
movaps xmm3, xmm1
shufps xmm1, xmm0, 129
shufps xmm2, xmm2, 0
movss xmm0, dword ptr [ecx+8]
mulps xmm4, xmm5
mulps xmm2, xmm1
shufps xmm0, xmm0, 0
addps xmm4, xmm2
mulps xmm0, xmm3
addps xmm4, xmm0
movss dword ptr [eax], xmm4
movhps qword ptr [eax+4], xmm4
}
StopRecordTime;
__asm ret
}

__declspec(naked) void PIII_Mult10_3x3_3x1(float *src1, float *src2, float *dst)
{
StartRecordTime;
__asm {
mov ecx, dword ptr [esp+8] ; src2
mov edx, dword ptr [esp+4] ; src1
mov eax, dword ptr [esp+0Ch] ; dst
movss xmm0, dword ptr [ecx]
movss xmm5, dword ptr [edx]
movhps xmm5, qword ptr [edx+4]
shufps xmm0, xmm0, 0
movss xmm1, dword ptr [ecx+4]
movss xmm3, dword ptr [edx+12]
movhps xmm3, qword ptr [edx+16]
shufps xmm1, xmm1, 0
mulps xmm0, xmm5
mulps xmm1, xmm3
movss xmm2, dword ptr [ecx+8]
shufps xmm2, xmm2, 0
movss xmm4, dword ptr [edx+24]
movhps xmm4, qword ptr [edx+28]
addps xmm0, xmm1
mulps xmm2, xmm4
addps xmm0, xmm2
movss dword ptr [eax], xmm0
movhps qword ptr [eax+4], xmm0
}
StopRecordTime;
__asm ret
}

__declspec(naked) void PIII_Mult00_4x4_4x1(float *src1, float *src2, float *dst)
{
StartRecordTime;
__asm {
mov ecx, dword ptr [esp+ 8] ; src2
mov edx, dword ptr [esp+ 4] ; src1
movlps xmm6, qword ptr [ecx ]
movlps xmm0, qword ptr [edx ]
shufps xmm6, xmm6, 0x44
movhps xmm0, qword ptr [edx+16]
mulps xmm0, xmm6
movlps xmm7, qword ptr [ecx+ 8]
movlps xmm2, qword ptr [edx+ 8]
shufps xmm7, xmm7, 0x44
movhps xmm2, qword ptr [edx+24]
mulps xmm2, xmm7
movlps xmm1, qword ptr [edx+32]
movhps xmm1, qword ptr [edx+48]
mulps xmm1, xmm6
movlps xmm3, qword ptr [edx+40]
addps xmm0, xmm2
movhps xmm3, qword ptr [edx+56]
mov eax, dword ptr [esp+12] ; dst
mulps xmm3, xmm7
movaps xmm4, xmm0
addps xmm1, xmm3
shufps xmm4, xmm1, 0x88
shufps xmm0, xmm1, 0xDD
addps xmm0, xmm4
movaps xmmword ptr [eax], xmm0
}
StopRecordTime;
__asm ret
}

__declspec(naked)
void PIII_Mult00_6x6_6x1(float *src1, float *src2, float *dst)
{
StartRecordTime;
__asm {
mov ebx, dword ptr [esp+ 4] ; src1
mov ecx, dword ptr [esp+ 8] ; src2
movlps xmm7, qword ptr [ecx]
movlps xmm6, qword ptr [ecx+8]
shufps xmm7, xmm7, 0x44
shufps xmm6, xmm6, 0x44
movlps xmm0, qword ptr [ebx ]
movhps xmm0, qword ptr [ebx+ 24]
mulps xmm0, xmm7
movlps xmm3, qword ptr [ebx+ 8]
movhps xmm3, qword ptr [ebx+ 32]
mulps xmm3, xmm6
movlps xmm1, qword ptr [ebx+ 48]
movhps xmm1, qword ptr [ebx+ 72]
mulps xmm1, xmm7
movlps xmm2, qword ptr [ebx+ 96]
movhps xmm2, qword ptr [ebx+120]
mulps xmm2, xmm7
movlps xmm4, qword ptr [ebx+ 56]
movhps xmm4, qword ptr [ebx+ 80]
movlps xmm5, qword ptr [ebx+104]
movhps xmm5, qword ptr [ebx+128]
mulps xmm4, xmm6
movlps xmm7, qword ptr [ecx+16]
addps xmm0, xmm3
shufps xmm7, xmm7, 0x44
mulps xmm5, xmm6
addps xmm1, xmm4
movlps xmm3, qword ptr [ebx+ 16]
movhps xmm3, qword ptr [ebx+ 40]
addps xmm2, xmm5
movlps xmm4, qword ptr [ebx+ 64]
movhps xmm4, qword ptr [ebx+ 88]
mulps xmm3, xmm7
movlps xmm5, qword ptr [ebx+112]
movhps xmm5, qword ptr [ebx+136]
addps xmm0, xmm3
mulps xmm4, xmm7
mulps xmm5, xmm7
addps xmm1, xmm4
addps xmm2, xmm5
movaps xmm6, xmm0
shufps xmm0, xmm1, 0x88
shufps xmm6, xmm1, 0xDD
movaps xmm7, xmm2
shufps xmm7, xmm2, 0x88
mov eax, dword ptr [esp+12] ; dst
shufps xmm2, xmm2, 0xDD
addps xmm0, xmm6
addps xmm2, xmm7
movaps xmmword ptr [eax], xmm0
movlps qword ptr [eax+16], xmm2
}
StopRecordTime;
__asm ret
}

__declspec(naked) void PIII_Mult10_4x4_4x1(float *src1, float *src2, float *dst)
{
StartRecordTime;
__asm {
mov ecx, dword ptr [esp+8] ; src2
mov edx, dword ptr [esp+4] ; src1
movss xmm0, dword ptr [ecx]
mov eax, dword ptr [esp+0Ch] ; dst
shufps xmm0, xmm0, 0
movss xmm1, dword ptr [ecx+4]
mulps xmm0, xmmword ptr [edx]
shufps xmm1, xmm1, 0
movss xmm2, dword ptr [ecx+8]
mulps xmm1, xmmword ptr [edx+16]
shufps xmm2, xmm2, 0
movss xmm3, dword ptr [ecx+12]
mulps xmm2, xmmword ptr [edx+32]
shufps xmm3, xmm3, 0
addps xmm0, xmm1
mulps xmm3, xmmword ptr [edx+48]
addps xmm2, xmm3
addps xmm0, xmm2
movaps xmmword ptr [eax], xmm0
}
StopRecordTime;
__asm ret
}

int Ra;
int Ca;
int Rb;
int Cb;
int StrideA; // Stride form one row of A to the next (in bytes)
int StrideB; // Stride form one row of B to the next (in bytes)

void MatrixMult(float *MatrixA, float *MatrixB, float *MatrixO)
{
StartRecordTime;
__asm {
pushad
Matrix_of_Results_Setup:
mov ecx, 0 ; Counter for rows in A - Ra
Row_of_Results_Loop:
mov ebx, 0 ; Counter for columns in B - Cb
Dot_Product_Setup:
mov eax, 0 ; Counter for single dot product - Ca or Rb
mov esi, MatrixA ; Load pointer to An0
mov edi, MatrixB ; Load pointer to B00
lea edi, [edi+ebx*4] ; Adjust pointer horizontally to correct batch of 24
xorps xmm2, xmm2 ; zero out accumulators for pass of 24 results
xorps xmm3, xmm3
xorps xmm4, xmm4
xorps xmm5, xmm5
xorps xmm6, xmm6
xorps xmm7, xmm7
Dot_Product_Loop:
mov edx, [esi+eax*4]
shl edx, 1
cmp edx, 0
je Sparse_Entry_Escape
movss xmm0, [esi+eax*4]
shufps xmm0, xmm0, 0x0
movaps xmm1, [edi]
mulps xmm1, xmm0
addps xmm2, xmm1
movaps xmm1, [edi+16]
mulps xmm1, xmm0
addps xmm3, xmm1
movaps xmm1, [edi+32]
mulps xmm1, xmm0
addps xmm4, xmm1
movaps xmm1, [edi+48]
mulps xmm1, xmm0
addps xmm5, xmm1
movaps xmm1, [edi+64]
mulps xmm1, xmm0
addps xmm6, xmm1
movaps xmm1, [edi+80]
mulps xmm1, xmm0
addps xmm7, xmm1
Sparse_Entry_Escape:
add edi, StrideB ; Move down a row in B
inc eax
cmp eax, Ca ; Can compare to Ca or Rb since they must be equal
jl Dot_Product_Loop
; End_Dot_Product_Loop
mov eax, MatrixO ; Load pointer to On0
lea eax, [eax+ebx*4] ; Adjust pointer horizontally to correct batch of 24
movaps [eax], xmm2 ; store to Output
movaps [eax+16], xmm3
movaps [eax+32], xmm4
movaps [eax+48], xmm5
movaps [eax+64], xmm6
movaps [eax+80], xmm7
add ebx, 24 ; Move over to next batch of 24
cmp ebx, Cb ; Check to see if row is complete
jl Dot_Product_Setup
; End_Row_of_Results_Loop
mov eax, MatrixA
add eax, StrideA
mov MatrixA, eax
mov eax, MatrixO
add eax, StrideB
mov MatrixO, eax
inc ecx
cmp ecx, Ra
jl Row_of_Results_Loop
; End_Matrix_Matrix_Multiply_Loop
popad
}
StopRecordTime;
}



char* testname;
_MM_ALIGN16 float m31[] = { 11,12,13,21,22,23,31,32,33};

_MM_ALIGN16 float m32[] = { 1,2,3,0,1,2,-1,0,1};

_MM_ALIGN16 float m33[] = {-2,34,70, -2,64,130, -2,94,190};
_MM_ALIGN16 float v4[] = {0, 1, 2, 3};
_MM_ALIGN16 float m41[] = { 11,12,13,14,21,22,23,24,31,32,33,34,41,42,43,44};
_MM_ALIGN16 float m42[] = { 1,2,3,4,0,1,2,3,-1,0,1,2,0,1,2,3};
_MM_ALIGN16 float m43[] = {-2,48,98,148,-2,88,178,268,-2,128,258,388,-2,168,338,508};
_MM_ALIGN16 float m44[16];
_MM_ALIGN16 float m45[16];
_MM_ALIGN16 float m46[16];
_MM_ALIGN16 float m61[] = { 11,12,13,14,15,16,21,22,23,24,25,26,31,32,33,34,35,36,41,42,43,44,45,46,51,52,53,54,55,56,61,62,63,64,65,66};
_MM_ALIGN16 float m62[] = { 1,2,3,4,5,6,
0,1,2,3,4,5,
-1,0,1,2,3,4,
0,1,2,3,4,5,
1,2,3,4,5,6,
2,3,4,5,6,7};

_MM_ALIGN16 float m63[36];
_MM_ALIGN16 float m64[36];
_MM_ALIGN16 float m65[36];
_MM_ALIGN16 float m66[36];

void testbase() {
SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_TIME_CRITICAL);
StartRecordTime;
StopRecordTime;
ticks[0] = SMF::MATH::end_ClockCount - SMF::MATH::start_ClockCount;
StartRecordTime;
StopRecordTime;
ticks[1] = SMF::MATH::end_ClockCount - SMF::MATH::start_ClockCount;
StartRecordTime;
StopRecordTime;
ticks[2] = SMF::MATH::end_ClockCount - SMF::MATH::start_ClockCount;
StartRecordTime;
StopRecordTime;
ticks[3] = SMF::MATH::end_ClockCount - SMF::MATH::start_ClockCount;
StartRecordTime;
StopRecordTime;
ticks[4] = SMF::MATH::end_ClockCount - SMF::MATH::start_ClockCount;
StartRecordTime;
StopRecordTime;
ticks[5] = SMF::MATH::end_ClockCount - SMF::MATH::start_ClockCount;
StartRecordTime;
StopRecordTime;
ticks[6] = SMF::MATH::end_ClockCount - SMF::MATH::start_ClockCount;
StartRecordTime;
StopRecordTime;
ticks[7] = SMF::MATH::end_ClockCount - SMF::MATH::start_ClockCount;
StartRecordTime;
StopRecordTime;
ticks[8] = SMF::MATH::end_ClockCount - SMF::MATH::start_ClockCount;
StartRecordTime;
StopRecordTime;
ticks[9] = SMF::MATH::end_ClockCount - SMF::MATH::start_ClockCount;
SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_NORMAL);

 report("%i %i %i %i %i %i %i %i %i %i ",
 (int)ticks[0], (int)ticks[1], (int)ticks[2], (int)ticks[3], (int)ticks[4],
 (int)ticks[5], (int)ticks[6], (int)ticks[7], (int)ticks[8], (int)ticks[9]);

base = Duration(10);
report("Duration for %s:\t%i", testname, base);

}

void test_3x3_3x1_PII() {

START_MEASUREMENTS;
PII_Mult00_3x3_3x1(m31, m32, v4);
END_MEASUREMENTS;

}

void test_3x3_3x1_PIII() {

START_MEASUREMENTS;
PIII_Mult00_3x3_3x1(m31, m32, v4);
END_MEASUREMENTS;

}

void test_3x3T_3x1_PIII() {

START_MEASUREMENTS;
PIII_Mult10_3x3_3x1(m31, m32, v4);
END_MEASUREMENTS;

}

void test_4x4_4x1_PII() {

START_MEASUREMENTS;
PII_Mult00_4x4_4x1(m41, m42, v4);
END_MEASUREMENTS;

}

void test_4x4_4x1_PIII() {

START_MEASUREMENTS;
PIII_Mult00_4x4_4x1(m41, m42, v4);
END_MEASUREMENTS;

}

void test_4x4T_4x1_PIII() {

START_MEASUREMENTS;
PIII_Mult10_4x4_4x1(m41, m42, v4);
END_MEASUREMENTS;

}

void test_3x3_3x3_PII() {

START_MEASUREMENTS;
PII_Mult_3x3_3x3(m31, m32, m33);
END_MEASUREMENTS;

}

void test_3x3_3x3_PIII() {

START_MEASUREMENTS;
PIII_Mult00_3x3_3x3(m31, m32, m33);
END_MEASUREMENTS;

}

void test_4x4_4x4_PII() {

START_MEASUREMENTS;
PII_Mult00_4x4_4x4(m41, m42, m43);
END_MEASUREMENTS;

}

void test_4x4_4x4_PIII() {

START_MEASUREMENTS;
PIII_Mult00_4x4_4x4(m41, m42, m43);
END_MEASUREMENTS;

}

void test_6x6_6x1_PII() {

START_MEASUREMENTS;
PII_Mult00_6x6_6x1(m61, m62, m63);
END_MEASUREMENTS;

}

void test_6x6_6x1_PIII() {

START_MEASUREMENTS;
PIII_Mult00_6x6_6x1(m61, m62, m63);
END_MEASUREMENTS;

}

void test_6x6_6x6_PII() {

START_MEASUREMENTS;
PII_Mult00_6x6_6x6(m61, m62, m63);
END_MEASUREMENTS;

}

void test_6x6_6x6_PIII() {

START_MEASUREMENTS;
PIII_Mult00_6x6_6x6(m61, m62, m63);
END_MEASUREMENTS;

}

void test_Inverse_4x4_PII() {

START_MEASUREMENTS;
StartRecordTime;
SMF::MATH::PII_Inverse_4x4(m44);
StopRecordTime;
END_MEASUREMENTS;

}

void test_Inverse_4x4_PIII() {

START_MEASUREMENTS;
StartRecordTime;
SMF::MATH::PIII_Inverse_4x4(m45);
StopRecordTime;
END_MEASUREMENTS;

}

void test_Inverse_6x6_PII() {

START_MEASUREMENTS;
StartRecordTime;
SMF::MATH::PII_Inverse_6x6(m64);
StopRecordTime;
END_MEASUREMENTS;

}

void test_InverseG_6x6_PIII() {

START_MEASUREMENTS;
StartRecordTime;
SMF::MATH::PIII_InverseG_6x6(m65);
StopRecordTime;
END_MEASUREMENTS;

}

void test_InverseS_6x6_PIII() {

START_MEASUREMENTS;
StartRecordTime;
SMF::MATH::PIII_InverseS_6x6(m66);
StopRecordTime;
END_MEASUREMENTS;

}

void testMult_4x4_PIII() {

Ra = 4;
Ca = 4;
Rb = 4;
Cb = 4;
StrideA = (((Ca+3)>>2)<<4);
StrideB = (((Cb+3)>>2)<<4);
START_MEASUREMENTS;
MatrixMult(m41, m42, m43);
END_MEASUREMENTS;

}

//This is for the normal PC
//The programm entry point
int main(int argc,char *argv[])
{

		//==============================
	SMF::Debug::setDebugAll();
	SMF::Debug::setFilename("debug_math.txt");
	SMF::Debug::setDebugMethod(SMF::Debug::File);

	//===============================
// We are looking for the best value among SAMPLES to

// eliminate cache delays and effects of cpuid variable timing.

testname = "rdtsc base";

testbase();

testname = "3x3 * 3x1 (PII)";

test_3x3_3x1_PII();

testname = "Transpose(3x3) * 3x1 (PIII)";

test_3x3T_3x1_PIII();

testname = "3x3 * 3x1 (PIII)";

test_3x3_3x1_PIII();

testname = "4x4 * 4x1 (PII)";

test_4x4_4x1_PII();

testname = "Transpose(4x4) * 4x1 (PIII)";

test_4x4T_4x1_PIII();

testname = "4x4 * 4x1 (PIII)";

test_4x4_4x1_PIII();

testname = "3x3 * 3x3 (PII)";

test_3x3_3x3_PII();

testname = "3x3 * 3x3 (PIII)";

test_3x3_3x3_PIII();

testname = "4x4 * 4x4 (PII)";

test_4x4_4x4_PII();

testname = "4x4 * 4x4 (PIII)";

test_4x4_4x4_PIII();

testname = "6x6 * 6x1 (PII)";

test_6x6_6x1_PII();

testname = "6x6 * 6x1 (PIII)";

test_6x6_6x1_PIII();

testname = "6x6 * 6x6 (PII)";

test_6x6_6x6_PII();

testname = "6x6 * 6x6 (PIII)";

test_6x6_6x6_PIII();

testname = "4x4 * 4x4 (general case, PIII)";

testMult_4x4_PIII();

int i;

for(i = 0 ; i < 16; i++)

m44[i] = m45[i] = (float)rand() / RAND_MAX;

for(i = 0 ; i < 36; i++)

m64[i] = m65[i] = m66[i] = (float)rand() / RAND_MAX;

testname = "Inverse 4x4 special (PII)";

test_Inverse_4x4_PII();

testname = "Inverse 4x4 (PIII)";

test_Inverse_4x4_PIII();

testname = "Inverse 6x6 special (PII)";

test_Inverse_6x6_PII();

testname = "Inverse 6x6 generic (PIII)";

test_InverseG_6x6_PIII();

testname = "Inverse 6x6 special (PIII)";

test_InverseS_6x6_PIII();

//scanf("%i");

// Test inverse.

#define zero(x) ( fabs(x) < 1e-4 ? 1:0 )

for( i=0; i<16; i++ )
 if( !zero( m44[i] - m45[i] ) )

 break;

 if( i<16 )

 report("Test PIII_Invert_4x4 fail.");

 else

 report("Test PIII_Invert_4x4 passed.");

//

 for( i=0; i<36; i++ )

 if( !zero( m64[i] - m65[i] ) )

 break;

 if( i<36 )

 report("Test PIII_InvertG_6x6 fail.");

 else

 report("Test PIII_InvertG_6x6 passed.");



 for( i=0; i<36; i++ )
 if( !zero( m64[i] - m66[i] ) )

break;

if( i<36 )

report("Test PIII_InvertS_6x6 fail.");

 else

 report("Test PIII_InvertS_6x6 passed.");

// Verify multiplications:

// 74, 134, 194 for 3x3 * 3x1

// 130, 230, 330, 430 for 4x4 * 4x1

 report("%4.1f %4.1f %4.1f %4.1f", v4[0], v4[1], v4[2], v4[3]);

 report("%4.1f %4.1f %4.1f", m33[0], m33[1], m33[2]); // -2,34,70

 report("%4.1f %4.1f %4.1f", m33[3], m33[4], m33[5]); // -2,64,130

 report("%4.1f %4.1f %4.1f", m33[6], m33[7], m33[8]); // -2,94,190

 report("%4.1f %4.1f %4.1f %4.1f", m43[0], m43[1], m43[2], m43[3]); // -2,48,98,148

 report("%4.1f %4.1f %4.1f %4.1f", m43[4], m43[5], m43[6], m43[7]); // -2,88,178,268

 report("%4.1f %4.1f %4.1f %4.1f", m43[8], m43[9], m43[10],m43[11]); // -2,128,258,388
 report("%4.1f %4.1f %4.1f %4.1f", m43[12],m43[13],m43[14],m43[15]); // -2,168,338,508
 report("%4.1f %4.1f %4.1f %4.1f %4.1f %4.1f",
 m63[0], m63[1], m63[2] ,m63[3] ,m63[4] ,m63[5]); // 45,126,207,288,369,450
 report("%4.1f %4.1f %4.1f %4.1f %4.1f %4.1f",
 m63[6], m63[7], m63[8] ,m63[9] ,m63[10],m63[11]); // 75,216,357,498,639,780
 report("%4.1f %4.1f %4.1f %4.1f %4.1f %4.1f",
 m63[12],m63[13],m63[14],m63[15],m63[16],m63[17]); // 105,306,507,708,909,1110
 report("%4.1f %4.1f %4.1f %4.1f %4.1f %4.1f",
 m63[18],m63[19],m63[20],m63[21],m63[22],m63[23]); // 135,396,657,918,1179,1440
 report("%4.1f %4.1f %4.1f %4.1f %4.1f %4.1f",
 m63[24],m63[25],m63[26],m63[27],m63[28],m63[29]); // 165,486,807,1128,1449,1770
 report("%4.1f %4.1f %4.1f %4.1f %4.1f %4.1f",
 m63[30],m63[31],m63[32],m63[33],m63[34],m63[35]); // 195,576,957,1338,1719,2100
	
    return 0;

}



