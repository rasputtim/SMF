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
#ifndef __SMF_BLOCK_ALOC
#define __SMF_BLOCK_ALOC

#include "../SMF_Config.h"

namespace SMF{
namespace Util{

/**
 * \class CBlockAlloc
 *
 * \ingroup SMF_Util
 * 
 * \if pt_br
 * \brief	Alocador de memória memória em blocos
 * \note	Todos os objetos "type" do template são construídos.
 *	Contudo, o construtor não é chamado para objetos reutilizados
 * \elseif us_en
 * \brief 	Block based allocator for fixed size objects.
 *
 *	\note	All objects of the 'type' are properly constructed.
 *	However, the constructor is not called for re-used objects.
 *
 * \endif
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
template<class type, int blockSize>
class SMF_API CBlockAlloc {
public:
							CBlockAlloc();
							~CBlockAlloc();

	void					shutdown();

	type *					alloc();
	void					free( type *element );

	int						getTotalCount() const { return total; }
	int						getAllocCount() const { return active; }
	int						getFreeCount() const { return total - active; }

private:
	/**
	\brief estrutura que armazena um elemento do tipo 'type' e um ponteiro para o próximo elemento
	**/
	typedef struct element_s {
		struct element_s *	next;
		type				t;
	} element_t;
	/**
	\brief estrutura que armazena um bloco do tipo element_t e um ponteiro para o próximo bloco
	**/
	typedef struct block_s {
		element_t			elements[blockSize];
		struct block_s *	next;
	} block_t;

	block_t *				blocks;
	element_t *				freebl;
	int						total;
	int						active;
};

template<class type, int blockSize>
CBlockAlloc<type,blockSize>::CBlockAlloc() {
	blocks = NULL;
	freebl = NULL;
	total = active = 0;
}

template<class type, int blockSize>
CBlockAlloc<type,blockSize>::~CBlockAlloc() {
	shutdown();
}

template<class type, int blockSize>
type *CBlockAlloc<type,blockSize>::alloc() {
	if ( !freebl ) {
		block_t *block = new block_t;
		block->next = blocks;
		blocks = block;
		for ( int i = 0; i < blockSize; i++ ) {
			block->elements[i].next = freebl;
			freebl = &block->elements[i];
		}
		total += blockSize;
	}
	active++;
	element_t *element = freebl;
	freebl = freebl->next;
	element->next = NULL;
	return &element->t;
}

template<class type, int blockSize>
void CBlockAlloc<type,blockSize>::free( type *t ) {
	element_t *element = (element_t *)( ( (unsigned char *) t ) - ( (int) &((element_t *)0)->t ) );
	element->next = freebl;
	freebl = element;
	active--;
}

template<class type, int blockSize>
void CBlockAlloc<type,blockSize>::shutdown() {
	while( blocks ) {
		block_t *block = blocks;
		blocks = blocks->next;
		delete block;
	}
	blocks = NULL;
	freebl = NULL;
	total = active = 0;
}

} //end Util
} //end SMF
#endif // __SMF_BLOCK_ALOC