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

#ifndef _SMF__MATH_RANDOM_H__
#define _SMF__MATH_RANDOM_H__
#include "../SMF_Config.h"
#include "../math/SMF_Math.h"

namespace SMF {
using namespace MATH;
namespace Util{

/**
 * \class 	CRandom
 *
 * \ingroup SMF_Util
 *
 * \if pt_br
 * \brief	Gerador ramdômico de números
 *
 * \elseif us_en
 * \brief 	random number generator
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
class SMF_API CRandom {
public:
						CRandom( int seed = 0 );

	void				setSeed( int seed );
	int					GetSeed() const;
	/// random integer in the range [0, MAX_RAND]
	int					randomInt();		
	/// random integer in the range [0, max[
	int					randomInt( int max );		
	/// random number in the range [0.0f, 1.0f]
	float				randomFloat();		
	// random number in the range [-1.0f, 1.0f]
	float				randomFloatUnsigned();		

	static const int	MAX_RAND = 0x7fff;

private:
	int					seed;
};

SMF_INLINE_FORCED CRandom::CRandom( int seed ) {
	this->seed = seed;
}

SMF_INLINE_FORCED void CRandom::setSeed( int seed ) {
	this->seed = seed;
}

SMF_INLINE_FORCED int CRandom::GetSeed() const {
	return seed;
}

SMF_INLINE_FORCED int CRandom::randomInt() {
	seed = 69069 * seed + 1;
	return ( seed & CRandom::MAX_RAND );
}

SMF_INLINE_FORCED int CRandom::randomInt( int max ) {
	if ( max == 0 ) {
		return 0;			// avoid divide by zero error
	}
	return randomInt() % max;
}

SMF_INLINE_FORCED float CRandom::randomFloat() {
	return ( randomInt() / ( float )( CRandom::MAX_RAND + 1 ) );
}

SMF_INLINE_FORCED float CRandom::randomFloatUnsigned() {
	return ( 2.0f * ( randomFloat() - 0.5f ) );
}



/**
 * \class 	CRandom2
 *
 * \ingroup SMF_Util
 *
 * \if pt_br
 * \brief	Gerador ramdômico de números
 *
 * \elseif us_en
 * \brief 	random number generator
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

class CRandom2 {
public:
							CRandom2( unsigned long seed = 0 );

	void					setSeed( unsigned long seed );
	unsigned long			GetSeed() const;
	/// random integer in the range [0, MAX_RAND]
	int						randomInt();
	/// random integer in the range [0, max]
	int						randomInt( int max );		
	/// random number in the range [0.0f, 1.0f]
	float					randomFloat();
	/// random number in the range [-1.0f, 1.0f]
	float					randomFloatUnsigned();		

	static const int		MAX_RAND = 0x7fff;

private:
	unsigned long			seed;

	static const unsigned long	IEEE_ONE = 0x3f800000;
	static const unsigned long	IEEE_MASK = 0x007fffff;
};

SMF_INLINE_FORCED CRandom2::CRandom2( unsigned long seed ) {
	this->seed = seed;
}

SMF_INLINE_FORCED void CRandom2::setSeed( unsigned long seed ) {
	this->seed = seed;
}

SMF_INLINE_FORCED unsigned long CRandom2::GetSeed() const {
	return seed;
}

SMF_INLINE_FORCED int CRandom2::randomInt() {
	seed = 1664525L * seed + 1013904223L;
	return ( (int) seed & CRandom2::MAX_RAND );
}

SMF_INLINE_FORCED int CRandom2::randomInt( int max ) {
	if ( max == 0 ) {
		return 0;		// avoid divide by zero error
	}
	return ( randomInt() >> ( 16 - CMath::bitsForInteger( max ) ) ) % max;
}

SMF_INLINE_FORCED float CRandom2::randomFloat() {
	unsigned long i;
	seed = 1664525L * seed + 1013904223L;
	i = CRandom2::IEEE_ONE | ( seed & CRandom2::IEEE_MASK );
	return ( ( *(float *)&i ) - 1.0f );
}

SMF_INLINE_FORCED float CRandom2::randomFloatUnsigned() {
	unsigned long i;
	seed = 1664525L * seed + 1013904223L;
	i = CRandom2::IEEE_ONE | ( seed & CRandom2::IEEE_MASK );
	return ( 2.0f * ( *(float *)&i ) - 3.0f );
}

} //end MATH
} //end SMF
#endif /* !__MATH_RANDOM_H__ */
