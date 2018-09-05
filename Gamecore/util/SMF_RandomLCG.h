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

#ifndef _SMF__MATH_RANDOM_LCG_H__
#define _SMF__MATH_RANDOM_LCG_H__
#include "../SMF_Config.h"
#include "../math/SMF_Math.h"

namespace SMF {
using namespace MATH;
namespace Util{

/**
 * \class CRandomLCG
 *
 * \ingroup SMF_Util
 * \if pt_br
 * \brief Gerador de números randômicos com o método do D.H. Lehmer's Linear Congruential (1949)
 * \elseif us_en
 * \brief A linear congruential random number generator.

	Uses D.H. Lehmer's Linear Congruential Method (1949) for generating random numbers.
	Supports both Multiplicative Congruential Method (increment==0) and
	Mixed Congruential Method (increment!=0)
	It is perhaps the simplest and fastest method to generate pseudo-random numbers on
	a computer. Per default uses the values for Minimal Standard CRandomLCG.
	http://en.wikipedia.org/wiki/Linear_congruential_generator
	http://www.math.rutgers.edu/~greenfie/currentcourses/sem090/pdfstuff/jp.pdf

	Pros:
	<ul>
	    <li> Easy to implement.
	    <li> Fast.
	</ul>

	Cons:
	<ul>
	    <li> NOT safe for cryptography because of the easily calculatable sequential
	        correlation between successive calls. A case study:
	        http://www.cigital.com/papers/download/developer_gambling.php

	    <li> Tends to have less random low-order bits (compared to the high-order bits)
	         Thus, NEVER do something like this:

	           sf_u32 numBetween1And10 = 1 + LCGRand.getInt() % 10;

	         Instead, take into account EVERY bit of the generated number, like this:

	           sf_u32 numBetween1And10 = 1 + (int)(10.0 * (double)LCGRand.getInt()
	                                                      /(LCGRand.MAX()+1.0));
	         or simply
	
	           sf_u32 numBetween1And10 = LCGRand.Float(1.f, 10.f);
	</ul>
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
   */
class SMF_API CRandomLCG
{
public:
	/// Initializes the generator from the current system clock.
	CRandomLCG();
	/// Initializes the generator using a custom seed.
	CRandomLCG(sf_u32 seed, sf_u32 multiplier = 69621,
		sf_u32 increment = 0, sf_u32 modulus = 0x7FFFFFFF /* 2^31 - 1 */)
	{
		this->seed(seed, multiplier, increment, modulus);
	}

	/// Reinitializes the generator to the new settings.
	void seed(sf_u32 seed, sf_u32 multiplier = 69621, sf_u32 increment = 0, sf_u32 modulus = 0x7FFFFFFF);

	/// Returns an integer in the range [0, maxInt()]
	sf_u32 getInt();

	/// Returns the biggest number the generator can yield. (Which is always modulus-1)
	sf_u32 maxInt() const { return modulus - 1; }

	/**
	\brief Returns an integer in the range [0, 2^32-1].
	\note The configurable modulus and increment are not used by this function, but are always increment == 0, modulus=2^32.
	**/
	sf_u32 intFast();

	/** Returns an integer in the range [a, b]
	\param a Lower bound, inclusive.
	\param b Upper bound, inclusive.
	\return An integer in the range [a, b] 
	*/
	int getInt(int a, int b);

	/// Returns a float in the range [0, 1[.
	float getFloat();

	/**
	\brief Returns a float in the range [0, 1].
	 \note This is much slower than Float()! Prefer that function instead if possible.
	*/
	float getFloat01Incl();

	/** 
	\brief Returns a float in the range ]-1, 1[.
	\note This function has one more bit of randomness compared to Float(), but has a theoretical bias
	towards 0.0, since floating point has two representations for 0 (+0 and -0).
	**/
	float getFloatNeg1_1();

	/** Returns a float in the range [a, b[.
	\param a Lower bound, inclusive.
	\param b Upper bound, exclusive.
	\return A float in the range [a, b[ 
	*/
	float getFloat(float a, float b);

private:
	sf_u32 multiplier;
	sf_u32 increment;
	sf_u32 modulus;

	sf_u32 lastNumber;
};

#ifdef MATH_QT_INTEROP
Q_DECLARE_METATYPE(CRandomLCG)
Q_DECLARE_METATYPE(CRandomLCG*)
#endif


} //end MATH
} //end SMF
#endif