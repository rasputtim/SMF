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


#include "structures/SMF_HashTableFast.h"
#include "math/SMF_Math.h"

int SMF:: CHashTableFast::INVALID_INDEX[1] = { -1 };




/*
================
CHashTableFast::init
================
*/
void SMF:: CHashTableFast::init( const int initialHashSize, const int initialIndexSize ) {
	SMF_ASSERT( MATH::CMath::isPowerOfTwo( initialHashSize ) );

	m_iHashSize = initialHashSize;
	m_pHash = INVALID_INDEX;
	m_iIndexSize = initialIndexSize;
	m_pIndexChain = INVALID_INDEX;
	m_iGranularity = DEFAULT_HASH_GRANULARITY;
	m_iHashMask = m_iHashSize - 1;
	m_iLookupMask = 0;
}

/*
================
CHashTableFast::allocate
================
*/
void SMF:: CHashTableFast::allocate( const int newHashSize, const int newIndexSize ) {
	SMF_ASSERT( MATH::CMath::isPowerOfTwo( newHashSize ) );

	freeMem();
	m_iHashSize = newHashSize;
	m_pHash = new int[m_iHashSize];
	memset( m_pHash, 0xff, m_iHashSize * sizeof( m_pHash[0] ) );
	m_iIndexSize = newIndexSize;
	m_pIndexChain = new int[m_iIndexSize];
	memset( m_pIndexChain, 0xff, m_iIndexSize * sizeof( m_pIndexChain[0] ) );
	m_iHashMask = m_iHashSize - 1;
	m_iLookupMask = -1;
}

/*
================
CHashTableFast::freeMem
================
*/
void SMF:: CHashTableFast::freeMem() {
	if ( m_pHash != INVALID_INDEX ) {
		delete[] m_pHash;
		m_pHash = INVALID_INDEX;
	}
	if ( m_pIndexChain != INVALID_INDEX ) {
		delete[] m_pIndexChain;
		m_pIndexChain = INVALID_INDEX;
	}
	m_iLookupMask = 0;
}

/*
================
CHashTableFast::resizeIndex
================
*/
void SMF:: CHashTableFast::resizeIndex( const int newIndexSize ) {
	int *oldIndexChain, mod, newSize;

	if ( newIndexSize <= m_iIndexSize ) {
		return;
	}

	mod = newIndexSize % m_iGranularity;
	if ( !mod ) {
		newSize = newIndexSize;
	} else {
		newSize = newIndexSize + m_iGranularity - mod;
	}

	if ( m_pIndexChain == INVALID_INDEX ) {
		m_iIndexSize = newSize;
		return;
	}

	oldIndexChain = m_pIndexChain;
	m_pIndexChain = new int[newSize];
	memcpy( m_pIndexChain, oldIndexChain, m_iIndexSize * sizeof(int) );
	memset( m_pIndexChain + m_iIndexSize, 0xff, (newSize - m_iIndexSize) * sizeof(int) );
	delete[] oldIndexChain;
	m_iIndexSize = newSize;
}

/*
================
CHashTableFast::getSpread
================
*/
int SMF:: CHashTableFast::getSpread() const {
	int i, index, totalItems, *numHashItems, average, error, e;

	if ( m_pHash == INVALID_INDEX ) {
		return 100;
	}

	totalItems = 0;
	numHashItems = new int[m_iHashSize];
	for ( i = 0; i < m_iHashSize; i++ ) {
		numHashItems[i] = 0;
		for ( index = m_pHash[i]; index >= 0; index = m_pIndexChain[index] ) {
			numHashItems[i]++;
		}
		totalItems += numHashItems[i];
	}
	// if no items in m_pHash
	if ( totalItems <= 1 ) {
		delete[] numHashItems;
		return 100;
	}
	average = totalItems / m_iHashSize;
	error = 0;
	for ( i = 0; i < m_iHashSize; i++ ) {
		e = CMath::abs( numHashItems[i] - average );
		if ( e > 1 ) {
			error += e - 1;
		}
	}
	delete[] numHashItems;
	return 100 - (error * 100 / totalItems);
}

