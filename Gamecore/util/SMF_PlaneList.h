/*
  SMF -  Super Math Fabric  
  Copyright (C) 2014 Rasputtim <Rasputtim@hotmail.com>

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


#ifndef _SMF__PLANELIST_H__
#define _SMF__PLANELIST_H__
#include "../SMF_Config.h"
#include "../structures/SMF_List.h"
#include "../geometry/SMF_Plane.h"
#include "../structures/SMF_HashTableFast.h"
namespace SMF {
using namespace MATH;
namespace Util{
/**
 * \class CPlaneList
 *
 * \ingroup SMF_Data_Structures
 * \if pt_br
 * \brief Implementa uma listade Planos (CPlanes)
 * \elseif us_en
 * \brief CPlane List
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
   */
class SMF_API CPlaneList : public CList<CPlane> {
public:

	void					clear() { CList<CPlane>::clear(); m_hash.freeMem(); }

	int						findPlane( const CPlane &plane, const float normalEps, const float distEps );

private:
	CHashTableFast				m_hash;
};

SMF_INLINE_FORCED int CPlaneList::findPlane( const CPlane &plane, const float normalEps, const float distEps ) {
	int i, border, hashKey;

	SMF_ASSERT( distEps <= 0.125f );

	hashKey = (int)( CMath::fabs( plane.getDist() ) * 0.125f );
	for ( border = -1; border <= 1; border++ ) {
		for ( i = m_hash.getFirst( hashKey + border ); i >= 0; i = m_hash.getNext( i ) ) {
			if ( (*this)[i].compare( plane, normalEps, distEps ) ) {
				return i;
			}
		}
	}

	if ( plane.getType() >= PLANETYPE_NEGX && plane.getType() < PLANETYPE_TRUEAXIAL ) {
		append( -plane );
		m_hash.add( hashKey, getNum()-1 );
		append( plane );
		m_hash.add( hashKey, getNum()-1 );
		return ( getNum() - 1 );
	}
	else {
		append( plane );
		m_hash.add( hashKey, getNum()-1 );
		append( -plane );
		m_hash.add( hashKey, getNum()-1 );
		return ( getNum() - 2 );
	}
}

} //end MATH
} // end SMF
#endif /* !__SMF_PLANELIST_H__ */
