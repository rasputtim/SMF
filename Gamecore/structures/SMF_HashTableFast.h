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


#ifndef _SMF__HASHINDEX_H__
#define _SMF__HASHINDEX_H__
#include "../util/SMF_HashUtils.h"


#define DEFAULT_HASH_SIZE			1024
#define DEFAULT_HASH_GRANULARITY	1024
#include "../SMF_Config.h"
#include "../math/SMF_Matriz.h"

namespace SMF {
/**
 * \class CHashTableFast
 *
 * \ingroup SMF_Data_Structures
 * 
 * \if pt_br
 * \brief	Implementa uma Tabela de hash rápida
 * \note	Não aloca memória até que o primeiro par chave/indice seja adicionado
 * \elseif us_en
 * \brief 	Fast hash table for indexes and arrays.
 * \note	Does not allocate memory until the first key/index pair is added.

 * \endif
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \see http://en.wikipedia.org/wiki/Self-balancing_binary_search_tree
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CHashTableFast {
public:
					CHashTableFast();
					CHashTableFast( const int initialHashSize, const int initialIndexSize );
					~CHashTableFast();

					/// returns total size of allocated memory
	size_t			getAllocatedMem() const;
					/// returns total size of allocated memory including size of hash index type
	size_t			getSize() const;

	CHashTableFast &	operator=( const CHashTableFast &other );
					/// add an index to the hash, assumes the index has not yet been added to the hash
	void			add( const int key, const int index );
					/// remove an index from the hash
	void			remove( const int key, const int index );
					/// get the first index from the hash, returns -1 if empty hash entry
	int				getFirst( const int key ) const;
					/// get the next index from the hash, returns -1 if at the end of the hash chain
	int				getNext( const int index ) const;
					/// insert an entry into the index and add it to the hash, increasing all indexes >= index
	void			insertIndex( const int key, const int index );
					/// remove an entry from the index and remove it from the hash, decreasing all indexes >= index
	void			removeIndex( const int key, const int index );
					/// clear the hash
	void			clear();
					/// clear and resize
	void			clear( const int newHashSize, const int newIndexSize );
					/// free allocated memory
	void			freeMem();
					/// get size of hash table
	int				getHashSize() const;
					/// get size of the index
	int				getIndexSize() const;
					/// set m_iGranularity
	void			setGranularity( const int newGranularity );
					/// force resizing the index, current hash table stays intact
	void			resizeIndex( const int newIndexSize );
					/// returns number in the range [0-100] representing the spread over the hash table
	int				getSpread() const;
					/// returns a key for a string
	int				generate_Key( const char *string, bool caseSensitive = true ) const;
					/// returns a key for a vector
	int				generate_Key( const MATH::CVec3D &v ) const;
					/// returns a key for two integers
	int				generate_Key( const int n1, const int n2 ) const;

private:
	int				m_iHashSize;
	int *			m_pHash;
	int				m_iIndexSize;
	int *			m_pIndexChain;
	int				m_iGranularity;
	int				m_iHashMask;
	int				m_iLookupMask;

	static int		INVALID_INDEX[1];

	void			init( const int initialHashSize, const int initialIndexSize );
	void			allocate( const int newHashSize, const int newIndexSize );
};


SMF_INLINE_FORCED CHashTableFast::CHashTableFast() {
	init( DEFAULT_HASH_SIZE, DEFAULT_HASH_SIZE );
}


SMF_INLINE_FORCED CHashTableFast::CHashTableFast( const int initialHashSize, const int initialIndexSize ) {
	init( initialHashSize, initialIndexSize );
}


SMF_INLINE_FORCED CHashTableFast::~CHashTableFast() {
	freeMem();
}


SMF_INLINE_FORCED size_t CHashTableFast::getAllocatedMem() const {
	return m_iHashSize * sizeof( int ) + m_iIndexSize * sizeof( int );
}


SMF_INLINE_FORCED size_t CHashTableFast::getSize() const {
	return sizeof( *this ) + getAllocatedMem();
}


SMF_INLINE_FORCED CHashTableFast &CHashTableFast::operator=( const CHashTableFast &other ) {
	m_iGranularity = other.m_iGranularity;
	m_iHashMask = other.m_iHashMask;
	m_iLookupMask = other.m_iLookupMask;

	if ( other.m_iLookupMask == 0 ) {
		m_iHashSize = other.m_iHashSize;
		m_iIndexSize = other.m_iIndexSize;
		freeMem();
	}
	else {
		if ( other.m_iHashSize != m_iHashSize || m_pHash == INVALID_INDEX ) {
			if ( m_pHash != INVALID_INDEX ) {
				delete[] m_pHash;
			}
			m_iHashSize = other.m_iHashSize;
			m_pHash = new int[m_iHashSize];
		}
		if ( other.m_iIndexSize != m_iIndexSize || m_pIndexChain == INVALID_INDEX ) {
			if ( m_pIndexChain != INVALID_INDEX ) {
				delete[] m_pIndexChain;
			}
			m_iIndexSize = other.m_iIndexSize;
			m_pIndexChain = new int[m_iIndexSize];
		}
		memcpy( m_pHash, other.m_pHash, m_iHashSize * sizeof( m_pHash[0] ) );
		memcpy( m_pIndexChain, other.m_pIndexChain, m_iIndexSize * sizeof( m_pIndexChain[0] ) );
	}

	return *this;
}


SMF_INLINE_FORCED void CHashTableFast::add( const int key, const int index ) {
	int h;

	SMF_ASSERT( index >= 0 );
	if ( m_pHash == INVALID_INDEX ) {
		allocate( m_iHashSize, index >= m_iIndexSize ? index + 1 : m_iIndexSize );
	}
	else if ( index >= m_iIndexSize ) {
		resizeIndex( index + 1 );
	}
	h = key & m_iHashMask;
	m_pIndexChain[index] = m_pHash[h];
	m_pHash[h] = index;
}


SMF_INLINE_FORCED void CHashTableFast::remove( const int key, const int index ) {
	int k = key & m_iHashMask;

	if ( m_pHash == INVALID_INDEX ) {
		return;
	}
	if ( m_pHash[k] == index ) {
		m_pHash[k] = m_pIndexChain[index];
	}
	else {
		for ( int i = m_pHash[k]; i != -1; i = m_pIndexChain[i] ) {
			if ( m_pIndexChain[i] == index ) {
				m_pIndexChain[i] = m_pIndexChain[index];
				break;
			}
		}
	}
	m_pIndexChain[index] = -1;
}


SMF_INLINE_FORCED int CHashTableFast::getFirst( const int key ) const {
	return m_pHash[key & m_iHashMask & m_iLookupMask];
}


SMF_INLINE_FORCED int CHashTableFast::getNext( const int index ) const {
	SMF_ASSERT( index >= 0 && index < m_iIndexSize );
	return m_pIndexChain[index & m_iLookupMask];
}


SMF_INLINE_FORCED void CHashTableFast::insertIndex( const int key, const int index ) {
	int i, max;

	if ( m_pHash != INVALID_INDEX ) {
		max = index;
		for ( i = 0; i < m_iHashSize; i++ ) {
			if ( m_pHash[i] >= index ) {
				m_pHash[i]++;
				if ( m_pHash[i] > max ) {
					max = m_pHash[i];
				}
			}
		}
		for ( i = 0; i < m_iIndexSize; i++ ) {
			if ( m_pIndexChain[i] >= index ) {
				m_pIndexChain[i]++;
				if ( m_pIndexChain[i] > max ) {
					max = m_pIndexChain[i];
				}
			}
		}
		if ( max >= m_iIndexSize ) {
			resizeIndex( max + 1 );
		}
		for ( i = max; i > index; i-- ) {
			m_pIndexChain[i] = m_pIndexChain[i-1];
		}
		m_pIndexChain[index] = -1;
	}
	add( key, index );
}


SMF_INLINE_FORCED void CHashTableFast::removeIndex( const int key, const int index ) {
	int i, max;

	remove( key, index );
	if ( m_pHash != INVALID_INDEX ) {
		max = index;
		for ( i = 0; i < m_iHashSize; i++ ) {
			if ( m_pHash[i] >= index ) {
				if ( m_pHash[i] > max ) {
					max = m_pHash[i];
				}
				m_pHash[i]--;
			}
		}
		for ( i = 0; i < m_iIndexSize; i++ ) {
			if ( m_pIndexChain[i] >= index ) {
				if ( m_pIndexChain[i] > max ) {
					max = m_pIndexChain[i];
				}
				m_pIndexChain[i]--;
			}
		}
		for ( i = index; i < max; i++ ) {
			m_pIndexChain[i] = m_pIndexChain[i+1];
		}
		m_pIndexChain[max] = -1;
	}
}


SMF_INLINE_FORCED void CHashTableFast::clear() {
	// only clear the m_pHash table because clearing the m_pIndexChain is not really needed
	if ( m_pHash != INVALID_INDEX ) {
		memset( m_pHash, 0xff, m_iHashSize * sizeof( m_pHash[0] ) );
	}
}


SMF_INLINE_FORCED void CHashTableFast::clear( const int newHashSize, const int newIndexSize ) {
	freeMem();
	m_iHashSize = newHashSize;
	m_iIndexSize = newIndexSize;
}


SMF_INLINE_FORCED int CHashTableFast::getHashSize() const {
	return m_iHashSize;
}


SMF_INLINE_FORCED int CHashTableFast::getIndexSize() const {
	return m_iIndexSize;
}


SMF_INLINE_FORCED void CHashTableFast::setGranularity( const int newGranularity ) {
	SMF_ASSERT( newGranularity > 0 );
	m_iGranularity = newGranularity;
}


SMF_INLINE_FORCED int CHashTableFast::generate_Key( const char *string, bool caseSensitive ) const {
	if ( caseSensitive ) {
		return ( Util::getHash( string ) & m_iHashMask );
	} else {
		return ( Util::getIHash( string ) & m_iHashMask );
	}
}


SMF_INLINE_FORCED int CHashTableFast::generate_Key( const MATH::CVec3D &v ) const {
	return ( (((int) v[0]) + ((int) v[1]) + ((int) v[2])) & m_iHashMask );
}


SMF_INLINE_FORCED int CHashTableFast::generate_Key( const int n1, const int n2 ) const {
	return ( ( n1 + n2 ) & m_iHashMask );
}
} //end SMF
#endif /* !__HASHINDEX_H__ */
