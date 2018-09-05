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

#ifndef _SMF__HEAP_H__
#define _SMF__HEAP_H__

#include "../SMF_Config.h"
#include "../structures/SMF_BTree.h"
#include "SMF_HeapBlockAloc.h"
#include "../math/SMF_Math.h"

#ifndef _WIN32
#include <sys/mman.h>
#else
#include <Windows.h>
#endif

namespace SMF{
namespace Util{

/*
===============================================================================

	Memory Management

	This is a replacement for the compiler heap code (i.e. "C" malloc() and
	free() calls). On average 2.5-3.0 times faster than MSVC malloc()/free().
	Worst case performance is 1.65 times faster and best case > 70 times.

===============================================================================
*/






void		mem_Init();
/**
\brief retorna verdadero se o controle de memmória foi inicializado, falso caso contrário
**/
bool        mem_isInitialized();
void		mem_Shutdown();
void		mem_EnableLeakTest( const char *name );
void		mem_ClearFrameStats();
void		mem_GetFrameStats( memoryStats_t &allocs, memoryStats_t &frees );
void		mem_GetStats( memoryStats_t &stats );
void		mem_Dump_f( const class CCMDLineArgs &args );
void		mem_DumpCompressed_f( const class CCMDLineArgs &args );
void		mem_AllocDefragBlock();

/* *Allocates the given amount of memory at the given alignment.
\note
An "aligned" pointer by definition means that the numeric value of the pointer is evenly divisible by N (where N is the desired alignment).
To check this, cast the pointer to an integer of suitable size, take the modulus N, and check whether the result is zero.
 In code:
 bool is_aligned(void *p, int ALIGN)
{
    return (int)p % ALIGN == 0;
}n
 To check for alignment of a power of 2 you can use:

((unsigned long)p & (ALIGN - 1)) == 0

This is simply a faster version of (p % ALIGN) == 0.

(If ALIGN is a constant your compiler will probably automatically use the faster version above.)

\warning to use this function you do not need to initialize the CHeap. this method use standart libc malloc

**/
void * mem_alignedMalloc(size_t size, size_t alignment);

/// Frees memory allocated by mem_alignedMalloc.
void mem_alignedFree(void *ptr);




#ifndef SMF_DEBUG_MEM

void *		mem_Alloc( const int size );
void *		mem_ClearedAlloc( const int size );
void		mem_Free( void *ptr );
char *		mem_CopyString( const char *in );
void *		mem_Alloc16( const int size );
void		mem_Free16( void *ptr );

#ifdef SMF_REDIRECT_NEWDELETE

SMF_INLINE void *operator new( size_t s ) {
	return mem_Alloc( s );
}
SMF_INLINE void operator delete( void *p ) {
	mem_Free( p );
}
SMF_INLINE void *operator new[]( size_t s ) {
	return mem_Alloc( s );
}
SMF_INLINE void operator delete[]( void *p ) {
	mem_Free( p );
}

#endif

#else /* SMF_DEBUG_MEM */

void *		mem_Alloc( const int size, const char *fileName, const int lineNumber );
void *		mem_ClearedAlloc( const int size, const char *fileName, const int lineNumber );
void		mem_Free( void *ptr, const char *fileName, const int lineNumber );
char *		mem_CopyString( const char *in, const char *fileName, const int lineNumber );
void *		mem_Alloc16( const int size, const char *fileName, const int lineNumber );
void		mem_Free16( void *ptr, const char *fileName, const int lineNumber );

#ifdef SMF_REDIRECT_NEWDELETE

SMF_INLINE void *operator new( size_t s, int t1, int t2, char *fileName, int lineNumber ) {
	return mem_Alloc( s, fileName, lineNumber );
}
SMF_INLINE void operator delete( void *p, int t1, int t2, char *fileName, int lineNumber ) {
	mem_Free( p, fileName, lineNumber );
}
SMF_INLINE void *operator new[]( size_t s, int t1, int t2, char *fileName, int lineNumber ) {
	return mem_Alloc( s, fileName, lineNumber );
}
SMF_INLINE void operator delete[]( void *p, int t1, int t2, char *fileName, int lineNumber ) {
	mem_Free( p, fileName, lineNumber );
}
SMF_INLINE void *operator new( size_t s ) {
	return mem_Alloc( s, "", 0 );
}
SMF_INLINE void operator delete( void *p ) {
	mem_Free( p, "", 0 );
}
SMF_INLINE void *operator new[]( size_t s ) {
	return mem_Alloc( s, "", 0 );
}
SMF_INLINE void operator delete[]( void *p ) {
	mem_Free( p, "", 0 );
}

#define ID_DEBUG_NEW						new( 0, 0, __FILE__, __LINE__ )
#undef new
#define new									ID_DEBUG_NEW

#endif

#define		mem_Alloc( size )				mem_Alloc( size, __FILE__, __LINE__ )
#define		mem_ClearedAlloc( size )		mem_ClearedAlloc( size, __FILE__, __LINE__ )
#define		mem_Free( ptr )					mem_Free( ptr, __FILE__, __LINE__ )
#define		mem_CopyString( s )				mem_CopyString( s, __FILE__, __LINE__ )
#define		mem_Alloc16( size )				mem_Alloc16( size, __FILE__, __LINE__ )
#define		mem_Free16( ptr )				mem_Free16( ptr, __FILE__, __LINE__ )

#endif /* SMF_DEBUG_MEM */



/**
 * \class CDynamicAlloc
 *
 * \ingroup SMF_Util
 *
 * \if pt_br
 * \brief Alocador de memória dinâmico
 * \note	Nenhum construtor é chamado para o "type" do template
 * \note	Os blocos de memória alocados estão sempre alinhados em 16 bytes
 * \elseif us_en
 * \brief 	Dynamic Memory allocator, simple wrapper for normal allocations which can
 *	        be interchanged with CDynamicBlockAlloc.
 *
 *	\note	No constructor is called for the 'type'.
 *	\note	allocated blocks are always 16 bytes aligned.
 * \endif
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
template<class type, int baseBlockSize, int minBlockSize>
class SMF_API CDynamicAlloc {
public:
									CDynamicAlloc();
									~CDynamicAlloc();

	void							init();
	void							shutdown();
	void							setFixedBlocks( int numBlocks ) {}
	void							setLockMemory( bool lock ) {}
	void							FreeEmptyBaseBlocks() {}

	type *							alloc( const int num );
	type *							resize( type *ptr, const int num );
	void							free( type *ptr );
	const char *					checkMemory( const type *ptr ) const;

	int								getNumBaseBlocks() const { return 0; }
	int								getBaseBlockMemory() const { return 0; }
	int								getNumUsedBlocks() const { return numUsedBlocks; }
	int								getUsedBlockMemory() const { return usedBlockMemory; }
	int								getNumFreeBlocks() const { return 0; }
	int								getFreeBlockMemory() const { return 0; }
	int								getNumEmptyBaseBlocks() const { return 0; }

private:
	/// number of used blocks
	int								numUsedBlocks;
	/// total memory in used blocks
	int								usedBlockMemory;

	int								numAllocs;
	int								numResizes;
	int								numFrees;

	void							clear();
};

template<class type, int baseBlockSize, int minBlockSize>
CDynamicAlloc<type, baseBlockSize, minBlockSize>::CDynamicAlloc() {
	clear();
}

template<class type, int baseBlockSize, int minBlockSize>
CDynamicAlloc<type, baseBlockSize, minBlockSize>::~CDynamicAlloc() {
	shutdown();
}

template<class type, int baseBlockSize, int minBlockSize>
void CDynamicAlloc<type, baseBlockSize, minBlockSize>::init() {
}

template<class type, int baseBlockSize, int minBlockSize>
void CDynamicAlloc<type, baseBlockSize, minBlockSize>::shutdown() {
	clear();
}

template<class type, int baseBlockSize, int minBlockSize>
type *CDynamicAlloc<type, baseBlockSize, minBlockSize>::alloc( const int num ) {
	numAllocs++;
	if ( num <= 0 ) {
		return NULL;
	}
	numUsedBlocks++;
	usedBlockMemory += num * sizeof( type );
	return mem_Alloc16( num * sizeof( type ) );
}

template<class type, int baseBlockSize, int minBlockSize>
type *CDynamicAlloc<type, baseBlockSize, minBlockSize>::resize( type *ptr, const int num ) {

	numResizes++;

	if ( ptr == NULL ) {
		return alloc( num );
	}

	if ( num <= 0 ) {
		free( ptr );
		return NULL;
	}

	SMF_ASSERT( 0 );
	return ptr;
}

template<class type, int baseBlockSize, int minBlockSize>
void CDynamicAlloc<type, baseBlockSize, minBlockSize>::free( type *ptr ) {
	numFrees++;
	if ( ptr == NULL ) {
		return;
	}
	mem_Free16( ptr );
}

template<class type, int baseBlockSize, int minBlockSize>
const char *CDynamicAlloc<type, baseBlockSize, minBlockSize>::checkMemory( const type *ptr ) const {
	return NULL;
}

template<class type, int baseBlockSize, int minBlockSize>
void CDynamicAlloc<type, baseBlockSize, minBlockSize>::clear() {
	numUsedBlocks = 0;
	usedBlockMemory = 0;
	numAllocs = 0;
	numResizes = 0;
	numFrees = 0;
}






//#define DYNAMIC_BLOCK_ALLOC_CHECK

/**
 * \class CDynamicBlock
 *
 * \ingroup SMF_Util
 *
 * \if pt_br
 * \brief	Template de bloco de memória para o alocador dinâmico
 * \elseif us_en
 * \brief 	dynamic block of memory for the dinâmic allocator.
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
template<class type>
class CDynamicBlock {
public:
	type *							GetMemory() const { return (type *)( ( (sf_u8 *) this ) + sizeof( CDynamicBlock<type> ) ); }
	int								getSize() const { return std::abs( size ); }
	void							setSize( int s, bool isBaseBlock ) { size = isBaseBlock ? -s : s; }
	bool							IsBaseBlock() const { return ( size < 0 ); }

#ifdef DYNAMIC_BLOCK_ALLOC_CHECK
	int								id[3];
	void *							allocator;
#endif

	int								size;					// size in bytes of the block
	CDynamicBlock<type> *			prev;					// previous memory block
	CDynamicBlock<type> *			next;					// next memory block
	CBTreeNode<CDynamicBlock<type>,int> *node;			// node in the B-Tree with free blocks
};

/**
 * \class CDynamicBlockAlloc
 *
 * \ingroup SMF_Util
 *
 * \if pt_br
 * \brief	Template Alocador dinâmico Rápido de blocos de memória
 * \note	Nenhum construtor para "type" do template é chamado
 * \note	Os blocos de memória alocados são alinhados em 16bytes
 * \elseif us_en
 * \brief 	Fast dynamic block allocator.
 *
 *	\note	No constructor is called for the 'type'.
 *	\note	allocated blocks are always 16 byte aligned.
 *
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
template<class type, int baseBlockSize, int minBlockSize>
class CDynamicBlockAlloc {
public:
									CDynamicBlockAlloc();
									~CDynamicBlockAlloc();

	void							init();
	void							shutdown();
	void							setFixedBlocks( int numBlocks );
	void							setLockMemory( bool lock );
	void							FreeEmptyBaseBlocks();

	type *							alloc( const int num );
	type *							resize( type *ptr, const int num );
	void							free( type *ptr );
	const char *					checkMemory( const type *ptr ) const;

	int								getNumBaseBlocks() const { return numBaseBlocks; }
	int								getBaseBlockMemory() const { return baseBlockMemory; }
	int								getNumUsedBlocks() const { return numUsedBlocks; }
	int								getUsedBlockMemory() const { return usedBlockMemory; }
	int								getNumFreeBlocks() const { return numFreeBlocks; }
	int								getFreeBlockMemory() const { return freeBlockMemory; }
	int								getNumEmptyBaseBlocks() const;

private:
	/// first block in list in order of increasing address
	CDynamicBlock<type> *			firstBlock;
	/// last block in list in order of increasing address
	CDynamicBlock<type> *			lastBlock;
	/// B-Tree with free memory blocks
	CBtreeBal<CDynamicBlock<type>,int,4> freeTree;
	/// allow base block allocations
	bool							allowAllocs;
	/// lock memory so it cannot get swapped out
	bool							lockMemory;

#ifdef DYNAMIC_BLOCK_ALLOC_CHECK
	int								blockId[3];
#endif

	int								numBaseBlocks;			// number of base blocks
	int								baseBlockMemory;		// total memory in base blocks
	int								numUsedBlocks;			// number of used blocks
	int								usedBlockMemory;		// total memory in used blocks
	int								numFreeBlocks;			// number of free blocks
	int								freeBlockMemory;		// total memory in free blocks

	int								numAllocs;
	int								numResizes;
	int								numFrees;

	void							clear();
	CDynamicBlock<type> *			AllocInternal( const int num );
	CDynamicBlock<type> *			ResizeInternal( CDynamicBlock<type> *block, const int num );
	void							FreeInternal( CDynamicBlock<type> *block );
	void							LinkFreeInternal( CDynamicBlock<type> *block );
	void							UnlinkFreeInternal( CDynamicBlock<type> *block );
	void							checkMemory() const;
};

template<class type, int baseBlockSize, int minBlockSize>
CDynamicBlockAlloc<type, baseBlockSize, minBlockSize>::CDynamicBlockAlloc() {
	clear();
}

template<class type, int baseBlockSize, int minBlockSize>
CDynamicBlockAlloc<type, baseBlockSize, minBlockSize>::~CDynamicBlockAlloc() {
	shutdown();
}

template<class type, int baseBlockSize, int minBlockSize>
void CDynamicBlockAlloc<type, baseBlockSize, minBlockSize>::init() {
	freeTree.init();
}

template<class type, int baseBlockSize, int minBlockSize>
void CDynamicBlockAlloc<type, baseBlockSize, minBlockSize>::shutdown() {
	CDynamicBlock<type> *block;

	for ( block = firstBlock; block != NULL; block = block->next ) {
		if ( block->node == NULL ) {
			FreeInternal( block );
		}
	}

	for ( block = firstBlock; block != NULL; block = firstBlock ) {
		firstBlock = block->next;
		SMF_ASSERT( block->IsBaseBlock() );
		if ( lockMemory ) {
			Sys_UnlockMemory( block, block->getSize() + (int)sizeof( CDynamicBlock<type> ) );
		}
		mem_Free16( block );
	}

	freeTree.shutdown();

	clear();
}

template<class type, int baseBlockSize, int minBlockSize>
void CDynamicBlockAlloc<type, baseBlockSize, minBlockSize>::setFixedBlocks( int numBlocks ) {
	CDynamicBlock<type> *block;

	for ( int i = numBaseBlocks; i < numBlocks; i++ ) {
		block = ( CDynamicBlock<type> * ) mem_Alloc16( baseBlockSize );
		if ( lockMemory ) {
			Sys_LockMemory( block, baseBlockSize );
		}
#ifdef DYNAMIC_BLOCK_ALLOC_CHECK
		memcpy( block->id, blockId, sizeof( block->id ) );
		block->allocator = (void*)this;
#endif
		block->setSize( baseBlockSize - (int)sizeof( CDynamicBlock<type> ), true );
		block->next = NULL;
		block->prev = lastBlock;
		if ( lastBlock ) {
			lastBlock->next = block;
		} else {
			firstBlock = block;
		}
		lastBlock = block;
		block->node = NULL;

		FreeInternal( block );

		numBaseBlocks++;
		baseBlockMemory += baseBlockSize;
	}

	allowAllocs = false;
}

template<class type, int baseBlockSize, int minBlockSize>
void CDynamicBlockAlloc<type, baseBlockSize, minBlockSize>::setLockMemory( bool lock ) {
	lockMemory = lock;
}

template<class type, int baseBlockSize, int minBlockSize>
void CDynamicBlockAlloc<type, baseBlockSize, minBlockSize>::FreeEmptyBaseBlocks() {
	CDynamicBlock<type> *block, *next;

	for ( block = firstBlock; block != NULL; block = next ) {
		next = block->next;

		if ( block->IsBaseBlock() && block->node != NULL && ( next == NULL || next->IsBaseBlock() ) ) {
			UnlinkFreeInternal( block );
			if ( block->prev ) {
				block->prev->next = block->next;
			} else {
				firstBlock = block->next;
			}
			if ( block->next ) {
				block->next->prev = block->prev;
			} else {
				lastBlock = block->prev;
			}
			if ( lockMemory ) {
				Sys_UnlockMemory( block, block->getSize() + (int)sizeof( CDynamicBlock<type> ) );
			}
			numBaseBlocks--;
			baseBlockMemory -= block->getSize() + (int)sizeof( CDynamicBlock<type> );
			mem_Free16( block );
		}
	}

#ifdef DYNAMIC_BLOCK_ALLOC_CHECK
	checkMemory();
#endif
}

template<class type, int baseBlockSize, int minBlockSize>
int CDynamicBlockAlloc<type, baseBlockSize, minBlockSize>::getNumEmptyBaseBlocks() const {
	int numEmptyBaseBlocks;
	CDynamicBlock<type> *block;

	numEmptyBaseBlocks = 0;
	for ( block = firstBlock; block != NULL; block = block->next ) {
		if ( block->IsBaseBlock() && block->node != NULL && ( block->next == NULL || block->next->IsBaseBlock() ) ) {
			numEmptyBaseBlocks++;
		}
	}
	return numEmptyBaseBlocks;
}

template<class type, int baseBlockSize, int minBlockSize>
type *CDynamicBlockAlloc<type, baseBlockSize, minBlockSize>::alloc( const int num ) {
	CDynamicBlock<type> *block;

	numAllocs++;

	if ( num <= 0 ) {
		return NULL;
	}

	block = AllocInternal( num );
	if ( block == NULL ) {
		return NULL;
	}
	block = ResizeInternal( block, num );
	if ( block == NULL ) {
		return NULL;
	}

#ifdef DYNAMIC_BLOCK_ALLOC_CHECK
	checkMemory();
#endif

	numUsedBlocks++;
	usedBlockMemory += block->getSize();

	return block->GetMemory();
}

template<class type, int baseBlockSize, int minBlockSize>
type *CDynamicBlockAlloc<type, baseBlockSize, minBlockSize>::resize( type *ptr, const int num ) {

	numResizes++;

	if ( ptr == NULL ) {
		return alloc( num );
	}

	if ( num <= 0 ) {
		free( ptr );
		return NULL;
	}

	CDynamicBlock<type> *block = ( CDynamicBlock<type> * ) ( ( (sf_u8 *) ptr ) - (int)sizeof( CDynamicBlock<type> ) );

	usedBlockMemory -= block->getSize();

	block = ResizeInternal( block, num );
	if ( block == NULL ) {
		return NULL;
	}

#ifdef DYNAMIC_BLOCK_ALLOC_CHECK
	checkMemory();
#endif

	usedBlockMemory += block->getSize();

	return block->GetMemory();
}

template<class type, int baseBlockSize, int minBlockSize>
void CDynamicBlockAlloc<type, baseBlockSize, minBlockSize>::free( type *ptr ) {

	numFrees++;

	if ( ptr == NULL ) {
		return;
	}

	CDynamicBlock<type> *block = ( CDynamicBlock<type> * ) ( ( (sf_u8 *) ptr ) - (int)sizeof( CDynamicBlock<type> ) );

	numUsedBlocks--;
	usedBlockMemory -= block->getSize();

	FreeInternal( block );

#ifdef DYNAMIC_BLOCK_ALLOC_CHECK
	checkMemory();
#endif
}

template<class type, int baseBlockSize, int minBlockSize>
const char *CDynamicBlockAlloc<type, baseBlockSize, minBlockSize>::checkMemory( const type *ptr ) const {
	CDynamicBlock<type> *block;

	if ( ptr == NULL ) {
		return NULL;
	}

	block = ( CDynamicBlock<type> * ) ( ( (sf_u8 *) ptr ) - (int)sizeof( CDynamicBlock<type> ) );

	if ( block->node != NULL ) {
		return "memory has been freed";
	}

#ifdef DYNAMIC_BLOCK_ALLOC_CHECK
	if ( block->id[0] != 0x11111111 || block->id[1] != 0x22222222 || block->id[2] != 0x33333333 ) {
		return "memory has invalid id";
	}
	if ( block->allocator != (void*)this ) {
		return "memory was allocated with different allocator";
	}
#endif

	/* base blocks can be larger than baseBlockSize which can cause this code to fail
	CDynamicBlock<type> *base;
	for ( base = firstBlock; base != NULL; base = base->next ) {
		if ( base->IsBaseBlock() ) {
			if ( ((int)block) >= ((int)base) && ((int)block) < ((int)base) + baseBlockSize ) {
				break;
			}
		}
	}
	if ( base == NULL ) {
		return "no base block found for memory";
	}
	*/

	return NULL;
}

template<class type, int baseBlockSize, int minBlockSize>
void CDynamicBlockAlloc<type, baseBlockSize, minBlockSize>::clear() {
	firstBlock = lastBlock = NULL;
	allowAllocs = true;
	lockMemory = false;
	numBaseBlocks = 0;
	baseBlockMemory = 0;
	numUsedBlocks = 0;
	usedBlockMemory = 0;
	numFreeBlocks = 0;
	freeBlockMemory = 0;
	numAllocs = 0;
	numResizes = 0;
	numFrees = 0;

#ifdef DYNAMIC_BLOCK_ALLOC_CHECK
	blockId[0] = 0x11111111;
	blockId[1] = 0x22222222;
	blockId[2] = 0x33333333;
#endif
}

template<class type, int baseBlockSize, int minBlockSize>
CDynamicBlock<type> *CDynamicBlockAlloc<type, baseBlockSize, minBlockSize>::AllocInternal( const int num ) {
	CDynamicBlock<type> *block;
	int alignedBytes = ( num * sizeof( type ) + 15 ) & ~15;

	block = freeTree.findSmallestLargerEqual( alignedBytes );
	if ( block != NULL ) {
		UnlinkFreeInternal( block );
	} else if ( allowAllocs ) {
		int allocSize = MAX( baseBlockSize, alignedBytes + (int)sizeof( CDynamicBlock<type> ) );
		block = ( CDynamicBlock<type> * ) mem_Alloc16( allocSize );
		if ( lockMemory ) {
			Sys_LockMemory( block, baseBlockSize );
		}
#ifdef DYNAMIC_BLOCK_ALLOC_CHECK
		memcpy( block->id, blockId, sizeof( block->id ) );
		block->allocator = (void*)this;
#endif
		block->setSize( allocSize - (int)sizeof( CDynamicBlock<type> ), true );
		block->next = NULL;
		block->prev = lastBlock;
		if ( lastBlock ) {
			lastBlock->next = block;
		} else {
			firstBlock = block;
		}
		lastBlock = block;
		block->node = NULL;

		numBaseBlocks++;
		baseBlockMemory += allocSize;
	}

	return block;
}

template<class type, int baseBlockSize, int minBlockSize>
CDynamicBlock<type> *CDynamicBlockAlloc<type, baseBlockSize, minBlockSize>::ResizeInternal( CDynamicBlock<type> *block, const int num ) {
	int alignedBytes = ( num * sizeof( type ) + 15 ) & ~15;

#ifdef DYNAMIC_BLOCK_ALLOC_CHECK
	SMF_ASSERT( block->id[0] == 0x11111111 && block->id[1] == 0x22222222 && block->id[2] == 0x33333333 && block->allocator == (void*)this );
#endif

	// if the new size is larger
	if ( alignedBytes > block->getSize() ) {

		CDynamicBlock<type> *nextBlock = block->next;

		// try to annexate the next block if it's free
		if ( nextBlock && !nextBlock->IsBaseBlock() && nextBlock->node != NULL &&
				block->getSize() + (int)sizeof( CDynamicBlock<type> ) + nextBlock->getSize() >= alignedBytes ) {

			UnlinkFreeInternal( nextBlock );
			block->setSize( block->getSize() + (int)sizeof( CDynamicBlock<type> ) + nextBlock->getSize(), block->IsBaseBlock() );
			block->next = nextBlock->next;
			if ( nextBlock->next ) {
				nextBlock->next->prev = block;
			} else {
				lastBlock = block;
			}
		} else {
			// allocate a new block and copy
			CDynamicBlock<type> *oldBlock = block;
			block = AllocInternal( num );
			if ( block == NULL ) {
				return NULL;
			}
			memcpy( block->GetMemory(), oldBlock->GetMemory(), oldBlock->getSize() );
			FreeInternal( oldBlock );
		}
	}

	// if the unused space at the end of this block is large enough to hold a block with at least one element
	if (( block->getSize() - alignedBytes - (int)sizeof( CDynamicBlock<type> )) < (MAX( minBlockSize, (int)sizeof( type ) )) ) {
		return block;
	}

	CDynamicBlock<type> *newBlock;

	newBlock = ( CDynamicBlock<type> * ) ( ( (sf_u8 *) block ) + (int)sizeof( CDynamicBlock<type> ) + alignedBytes );
#ifdef DYNAMIC_BLOCK_ALLOC_CHECK
	memcpy( newBlock->id, blockId, sizeof( newBlock->id ) );
	newBlock->allocator = (void*)this;
#endif
	newBlock->setSize( block->getSize() - alignedBytes - (int)sizeof( CDynamicBlock<type> ), false );
	newBlock->next = block->next;
	newBlock->prev = block;
	if ( newBlock->next ) {
		newBlock->next->prev = newBlock;
	} else {
		lastBlock = newBlock;
	}
	newBlock->node = NULL;
	block->next = newBlock;
	block->setSize( alignedBytes, block->IsBaseBlock() );

	FreeInternal( newBlock );

	return block;
}

template<class type, int baseBlockSize, int minBlockSize>
void CDynamicBlockAlloc<type, baseBlockSize, minBlockSize>::FreeInternal( CDynamicBlock<type> *block ) {

	SMF_ASSERT( block->node == NULL );

#ifdef DYNAMIC_BLOCK_ALLOC_CHECK
	SMF_ASSERT( block->id[0] == 0x11111111 && block->id[1] == 0x22222222 && block->id[2] == 0x33333333 && block->allocator == (void*)this );
#endif

	// try to merge with a next free block
	CDynamicBlock<type> *nextBlock = block->next;
	if ( nextBlock && !nextBlock->IsBaseBlock() && nextBlock->node != NULL ) {
		UnlinkFreeInternal( nextBlock );
		block->setSize( block->getSize() + (int)sizeof( CDynamicBlock<type> ) + nextBlock->getSize(), block->IsBaseBlock() );
		block->next = nextBlock->next;
		if ( nextBlock->next ) {
			nextBlock->next->prev = block;
		} else {
			lastBlock = block;
		}
	}

	// try to merge with a previous free block
	CDynamicBlock<type> *prevBlock = block->prev;
	if ( prevBlock && !block->IsBaseBlock() && prevBlock->node != NULL ) {
		UnlinkFreeInternal( prevBlock );
		prevBlock->setSize( prevBlock->getSize() + (int)sizeof( CDynamicBlock<type> ) + block->getSize(), prevBlock->IsBaseBlock() );
		prevBlock->next = block->next;
		if ( block->next ) {
			block->next->prev = prevBlock;
		} else {
			lastBlock = prevBlock;
		}
		LinkFreeInternal( prevBlock );
	} else {
		LinkFreeInternal( block );
	}
}

template<class type, int baseBlockSize, int minBlockSize>
SMF_INLINE_FORCED void CDynamicBlockAlloc<type, baseBlockSize, minBlockSize>::LinkFreeInternal( CDynamicBlock<type> *block ) {
	block->node = freeTree.add( block, block->getSize() );
	numFreeBlocks++;
	freeBlockMemory += block->getSize();
}

template<class type, int baseBlockSize, int minBlockSize>
SMF_INLINE_FORCED void CDynamicBlockAlloc<type, baseBlockSize, minBlockSize>::UnlinkFreeInternal( CDynamicBlock<type> *block ) {
	freeTree.remove( block->node );
	block->node = NULL;
	numFreeBlocks--;
	freeBlockMemory -= block->getSize();
}

template<class type, int baseBlockSize, int minBlockSize>
void CDynamicBlockAlloc<type, baseBlockSize, minBlockSize>::checkMemory() const {
	CDynamicBlock<type> *block;

	for ( block = firstBlock; block != NULL; block = block->next ) {
		// make sure the block is properly linked
		if ( block->prev == NULL ) {
			SMF_ASSERT( firstBlock == block );
		} else {
			SMF_ASSERT( block->prev->next == block );
		}
		if ( block->next == NULL ) {
			SMF_ASSERT( lastBlock == block );
		} else {
			SMF_ASSERT( block->next->prev == block );
		}
	}
}
} //end Util
} //end SMF
#endif /* !__SMF_HEAP_H__ */
