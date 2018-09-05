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


#include "SMF_Config.h"
#include "util/SMF_RandomLCG.h"
#include "math/SMF_Math.h"
#include "exceptions/all.h"

namespace SMF {
using namespace MATH;
namespace Util{
	//todo: solve seed generation
CRandomLCG::CRandomLCG()
{
	//s seed(Clock::TickU32());
	seed(1349284);
}

/** If you want to give different parameters for the generator, you should remember that:
	- modulus should be prime
	- modulus, increment and the multiplier should all be relative primes (increment can be 0)
	- modulus should be greater than multiplier and increment

	Most often you can leave increment = 0.

	A list of widely used values:
	- Park and Miller (Minimal standard): mul = 16807 (7^5)             mod = 2^31 - 1 (2147483647 == 0x7FFFFFFF)
	- Park and Miller #2:                 mul = 48271                   mod = 2^31 - 1
	- Park and Miller #3:                 mul = 69621                   mod = 2^31 - 1
	- SIMSCRIPT:                          mul = 630360016               mod = 2^31 - 1
	- URN12:                              mul = 452807053               mod = 2^31

	Infamous examples (Don't use!):
	- Classical ANSI C                    mul = 1103515245  inc = 12345 mod = 2^31
	- RANDU                               mul = 65539                   mod = 2^31  */
void CRandomLCG::seed(sf_u32 seed, sf_u32 mul, sf_u32 inc, sf_u32 mod)
{
	//s SMF_ASSERT((seed != 0 || inc != 0) && "Initializing CRandomLCG with seed=0 && inc=0 results in an infinite series of 0s!");

	lastNumber = seed;
	multiplier = mul;
	increment = inc;
	modulus = mod;
	SMF_ASSERT(modulus != 0);
}

sf_u32 CRandomLCG::intFast()
{
// The configurable modulus and increment are not used by this function.
	sf_u32 mul = lastNumber * multiplier;
	lastNumber = mul - (mul <= lastNumber?1:0); // Whenever we overflow, flip by one to avoid even multiplier always producing even results, since modulus is even.
	return lastNumber;
}

sf_u32 CRandomLCG::getInt()
{
	SMF_ASSERT(modulus != 0);
	/// \todo Convert to using Schrage's method for approximate factorization. (Numerical Recipes in C)

	// Currently we cast everything to 64-bit to avoid overflow, which is quite dumb.

	// Create the new random number
//#ifdef WIN32
	sf_u64 newNum = ((sf_u64)lastNumber * (sf_u64)multiplier + (sf_u64)increment) % (sf_u64)modulus;
//	sf_u32 m = lastNumber * multiplier;
//	sf_u32 i = m + increment;
//	sf_u32 f = i & 0x7FFFFFFF;
//	sf_u32 m = (lastNumber * 214013 + 2531011) & 0x7FFFFFFF;
//	unsigned __int64 newNum = (lastNumber * multiplier + increment) & 0x7FFFFFFF;
//#else
	// On console platform, we rely on using smaller sequences.
//	unsigned long newNum = ((unsigned long)lastNumber * (unsigned long)multiplier + (unsigned long)increment) % (unsigned long)modulus;
//#endif
	// Save the newly generated random number to use as seed for the next one.
//	lastNumber = m;//(sf_u32)newNum;
	lastNumber = (sf_u32)newNum;
	return lastNumber;
}

int CRandomLCG::getInt(int a, int b)
{
	SMF_ASSERT(a <= b && "Error in range!");

//	return a + (int)(getInt() * MAX()/(b-a));
	int num = a + (int)(getFloat() * (b-a+1));
//	SMF_ASSERT(num >= a);
//	SMF_ASSERT(num <= b);
	///\todo Some bug here - the result is not necessarily in the proper range.
	if (num < a)
		num = a;
	if (num > b)
		num = b;
	return num;
}

float CRandomLCG::getFloat()
{
	sf_u32 i = ((sf_u32)getInt() & 0x007FFFFF /* random mantissa */) | 0x3F800000 /* fixed exponent */;
	float f = ReinterpretAsFloat(i); // f is now in range [1, 2[
	f -= 1.f; // Map to range [0, 1[
	return f;
}

float CRandomLCG::getFloat01Incl()
{
	for(int i = 0; i < 100; ++i)
	{
		sf_u32 val = (sf_u32)getInt() & 0x00FFFFFF;
		if (val > 0x800000)
			continue;
		else if (val == 0x800000)
			return 1.0f;
		else
		{
			val |= 0x3F800000;
			return ReinterpretAsFloat(val) - 1.f;
		}
	}
	return getFloat();
}

float CRandomLCG::getFloatNeg1_1()
{
	sf_u32 i = (sf_u32)getInt();
	sf_u32 one = ((i & 0x00800000) << 8) /* random sign bit */ | 0x3F800000 /* fixed exponent */;
	i = one | (i & 0x007FFFFF) /* random mantissa */;
	float f = ReinterpretAsFloat(i); // f is now in range ]-2, -1[ union [1, 2].
	float fone = ReinterpretAsFloat(one); // +/- 1, of same sign as f.
	return f - fone;
}

float CRandomLCG::getFloat(float a, float b)
{
	SMF_ASSERT(a <= b ); //&& "CRandomLCG::getFloat(a,b): Error in range: b < a!"

	return getFloat()*(b-a)+a;
}


} //end MATH
} //end SMF
