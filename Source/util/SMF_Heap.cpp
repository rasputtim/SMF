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

#include "util/SMF_Heap.h"
#include "util/SMF_Debug.h"
#include "math/SMF_Simd.h"
#include "exceptions/all.h"
namespace SMF {
namespace Util{



bool Sys_LockMemory(void *ptr, int bytes){
#ifdef _WIN32
  // http://msdn.microsoft.com/en-us/library/windows/desktop/aa366895%28v=vs.85%29.aspx
  /**Locks the specified regionf the process's virtual address space into physical memory, 
   * ensuring that subsequent access to the region will not incur a page fault.
   * **/
return (VirtualLock(ptr,(SIZE_T)bytes )!= FALSE );
#else
// http://www.informit.com/articles/article.aspx?p=23618&seqNum=9

return (mlock(ptr,(size_t)bytes) != -1);
#endif  
  
}

bool Sys_UnlockMemory( void *ptr, int bytes) {
#ifdef _WIN32
return (VirtualUnlock(ptr,(SIZE_T)bytes )!= FALSE );
#else
return (munlock(ptr,(size_t)bytes)!= -1);
#endif   
  
}



//===============================================================
//
//	CHeap
//
//===============================================================

#define SMALL_HEADER_SIZE		( (int) ( sizeof( sf_u8 ) + sizeof( sf_u8 ) ) )
#define MEDIUM_HEADER_SIZE		( (int) ( sizeof( mediumHeapEntry_s ) + sizeof( sf_u8 ) ) )
#define LARGE_HEADER_SIZE		( (int) ( sizeof( sf_u32 * ) + sizeof( sf_u8 ) ) )

#define ALIGN_SIZE( bytes )		( ( (bytes) + HEAP_ALIGN - 1 ) & ~(HEAP_ALIGN - 1) )
#define SMALL_ALIGN( bytes )	( ALIGN_SIZE( (bytes) + SMALL_HEADER_SIZE ) - SMALL_HEADER_SIZE )
#define MEDIUM_SMALLEST_SIZE	( ALIGN_SIZE( 256 ) + ALIGN_SIZE( MEDIUM_HEADER_SIZE ) )

/**
 * \class CHeap
 *
 * \ingroup SMF_Util
 * \if pt_br
 * \brief Implementa um Heap de memória 
 * \note	Classe utilizada apenas pelas funções de gerenciamento de memória.
            Você não pode utilizá-la diretamente
 * \elseif us_en
 * \brief	Memory Heap
 * \note Used by memory management functions. You can not use it directly
 * \endif
 * \see Util::mem_Init(), Util::mem_Shutdown()
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class CHeap {

public:
					CHeap();
					/// frees all associated data
					~CHeap();				
	/**
	\brief initializa the Heap
	**/
	void			init();					
	/**
	\brief alocate X bytes of memory
	\param bytes amount of memory to alocate
	**/
	void *			Allocate( const sf_u32 bytes );	
	/**
	/brief free memory
	\param p pointer to the memory to be freed
	**/
	void			free( void *p );				
	/**
	\brief alocate X bytes of memry with 16 bytes alignment
	\param bytes ammount of memory to be alocated
	**/
	void *			Allocate16( const sf_u32 bytes );
	/**
	\brief free memory alocated with 16 bytes aligned by Alocate16
	\param p pointer to the memory to be freed
	\see Allocate16
	**/
	void			Free16( void *p );				
	/**
	\briefreturns size of allocated memory block
	\param p	= pointer to memory block
	\note	size may not be the same as the size in the original
			allocation request (due to block alignment reasons).
	**/
	sf_u32			Msize( void *p );				 

	/**
	\brief dump contents of the heap
	**/
	void			Dump( void  );
	/**
	hack for huge renderbumps
	**/
	void 			AllocDefragBlock();	
	///return true if the CHeap is initialized
	static bool		isInitialized() {return initialized;}
private:
	/// memory alignment in bytes
	enum {
		HEAP_ALIGN = 8									
	};

	enum {
		///invalid allocation
		INVALID_ALLOC	= 0xdd,
		/// small allocation
		SMALL_ALLOC		= 0xaa,
		/// medium allocaction						
		MEDIUM_ALLOC	= 0xbb,
		/// large allocaction						
		LARGE_ALLOC		= 0xcc						
	};

    /// allocation page
	struct page_s {		
		/// data pointer to allocated memory
		void *				data;
		/// number of bytes of memory 'data' points to					
		sf_u32			dataSize;	
		/// next free page in same page manager				
		page_s *			next;					
		/// used only when allocated
		page_s *			prev;					
		/// this data used by the medium-size heap manager
		sf_u32			largestFree;				
		/// pointer to first free entry
		void *				firstFree;				
	};
	/// header fro the block with medium size of memory alocation
	struct mediumHeapEntry_s {
		/// pointer to page
		page_s *			page;					
		/// size of block
		sf_u32			size;						
		/// previous block
		mediumHeapEntry_s *	prev;					
		/// next block
		mediumHeapEntry_s *	next;					
		/// previous free block
		mediumHeapEntry_s *	prevFree;				
		/// next free block
		mediumHeapEntry_s *	nextFree;				
		/// non-zero if free block
		sf_u32			freeBlock;					
	};

	// variables
	void *			smallFirstFree[256/HEAP_ALIGN+1];	//! small heap allocator lists (for allocs of 1-255 bytes)
	page_s *		smallCurPage;					//! current page for small allocations
	sf_u32			smallCurPageOffset;				//! bytes offset in current page
	page_s *		smallFirstUsedPage;				//! first used page of the small heap manager

	page_s *		mediumFirstFreePage;			//! first partially free page
	page_s *		mediumLastFreePage;				//! last partially free page
	page_s *		mediumFirstUsedPage;			//! completely used page

	page_s *		largeFirstUsedPage;				//! first page used by the large heap manager

	page_s *		swapPage;

	sf_u32			pagesAllocated;					//! number of pages currently allocated
	sf_u32			pageSize;						//! size of one alloc page in bytes

	sf_u32			pageRequests;					//! page requests
	sf_u32			OSAllocs;						//! number of allocs made to the OS
	static bool				initialized;
	int				c_heapAllocRunningCount;

	void			*defragBlock;					//! a single huge block that can be allocated
													//! at startup, then freed when needed

	/**
	\brief allocates memory from the OS
	\param bytes	= page size in bytes
	\return pointer to page
	**/
	page_s *		AllocatePage( sf_u32 bytes );
	/**
	\brief frees a page back to the operating system
	\param p	= pointer to page
	**/
	void			FreePage( CHeap::page_s *p );
//===============================================================
//
//	small heap code
//
//===============================================================

	/**
	\brief allocate memory (1-255 bytes) from the small heap manager
	\param bytes = number of bytes to allocate
	\return returns pointer to allocated memory
	**/
	void *			SmallAllocate( sf_u32 bytes );	
	/*
	\brief frees a block of memory allocated by SmallAllocate() call
	\param data = pointer to block of memory
	**/
	void			SmallFree( void *ptr );
//===============================================================
//
//	medium heap code
//
//	Medium-heap allocated pages not returned to OS until heap destructor
//	called (re-used instead on subsequent medium-size malloc requests).
//
//===============================================================

/**
	\brief performs allocation using the medium heap manager from a given page
	\param p		= page
	\param  sizeNeeded	= number of bytes needed
	\return returns pointer to allocated memory

*/
	void *			MediumAllocateFromPage( CHeap::page_s *p, sf_u32 sizeNeeded );
	/**
	\brief allocate memory (256-32768 bytes) from medium heap manager
	\param bytes	= number of bytes to allocate
	\return pointer to allocated memory
	**/
	void *			MediumAllocate( sf_u32 bytes );	
	/**
	\brief frees a block allocated by the medium heap manager
	\param ptr	= pointer to data block
	**/
	void			MediumFree( void *ptr );		
//===============================================================
//
//	large heap code
//
//===============================================================

	/*
	\brief allocates large block of memory from the operating system directly
	\param bytes	= number of bytes to allocate
	\return returns pointer to allocated memory
	*/
	void *			LargeAllocate( sf_u32 bytes );	
	/**
	\brief frees a block of memory allocated by the 'large memory allocator'
	\param  p	= pointer to allocated memory
	**/
	void			LargeFree( void *ptr );
	/**
	\brief releases the swap page to OS
	**/
	void			ReleaseSwappedPages();
	/**
	\brief frees page to be used by the OS
	\param p	= page to free
	**/
	void			FreePageReal( CHeap::page_s *p );
};

bool	CHeap::initialized=false;



void CHeap::init () {
	OSAllocs			= 0;
	pageRequests		= 0;
	pageSize			= 65536 - sizeof( CHeap::page_s );
	pagesAllocated		= 0;								// reset page allocation counter

	largeFirstUsedPage	= NULL;								// init large heap manager
	swapPage			= NULL;

	memset( smallFirstFree, 0, sizeof(smallFirstFree) );	// init small heap manager
	smallFirstUsedPage	= NULL;
	smallCurPage		= AllocatePage( pageSize );
	SMF_ASSERT( smallCurPage );
	smallCurPageOffset	= SMALL_ALIGN( 0 );

	defragBlock = NULL;

	mediumFirstFreePage	= NULL;								// init medium heap manager
	mediumLastFreePage	= NULL;
	mediumFirstUsedPage	= NULL;

	c_heapAllocRunningCount = 0;
	initialized=true;
}


CHeap::CHeap() {
	init();
}


CHeap::~CHeap() {

	CHeap::page_s	*p;

	if ( smallCurPage ) {
		FreePage( smallCurPage );			/// free small-heap current allocation page
	}
	p = smallFirstUsedPage;					/// free small-heap allocated pages 
	while( p ) {
		CHeap::page_s *next = p->next;
		FreePage( p );
		p= next;
	}

	p = largeFirstUsedPage;					/// free large-heap allocated pages
	while( p ) {
		CHeap::page_s *next = p->next;
		FreePage( p );
		p = next;
	}

	p = mediumFirstFreePage;				/// free medium-heap allocated pages
	while( p ) {
		CHeap::page_s *next = p->next;
		FreePage( p );
		p = next;
	}

	p = mediumFirstUsedPage;				/// free medium-heap allocated completely used pages
	while( p ) {
		CHeap::page_s *next = p->next;
		FreePage( p );
		p = next;
	}

	ReleaseSwappedPages();			

	if ( defragBlock ) {
		free( defragBlock );
	}

	SMF_ASSERT( pagesAllocated == 0 );
}


void CHeap::AllocDefragBlock() {
	int		size = 0x40000000;

	if ( defragBlock ) {
		return;
	}
	while( 1 ) {
		defragBlock = malloc( size );
		if ( defragBlock ) {
			break;
		}
		size >>= 1;
	}
	Debug::debug(Debug::heap,__FUNCTION__) << "allocated a %i mb defrag block size: " << size / (1024*1024) << endl;
}


void *CHeap::Allocate( const sf_u32 bytes ) {
	if ( !bytes ) {
		return NULL;
	}
	c_heapAllocRunningCount++;

#if USE_LIBC_MALLOC
	return malloc( bytes );
#else
	if ( !(bytes & ~255) ) {
		return SmallAllocate( bytes );
	}
	if ( !(bytes & ~32767) ) {
		return MediumAllocate( bytes );
	}
	return LargeAllocate( bytes );
#endif
}


void CHeap::free( void *p ) {
	if ( !p ) {
		return;
	}
	c_heapAllocRunningCount--;

#if USE_LIBC_MALLOC
	free( p );
#else
	switch( ((sf_u8 *)(p))[-1] ) {
		case SMALL_ALLOC: {
			SmallFree( p );
			break;
		}
		case MEDIUM_ALLOC: {
			MediumFree( p );
			break;
		}
		case LARGE_ALLOC: {
			LargeFree( p );
			break;
		}
		default: {
			break;
		}
	}
#endif
}


void *CHeap::Allocate16( const sf_u32 bytes ) {
	sf_u8 *ptr, *alignedPtr;

	ptr = (sf_u8 *) malloc( bytes + 16 + 4 );
	if ( !ptr ) {
		if ( defragBlock ) {
			Debug::debug(Debug::heap,__FUNCTION__) << "Freeing defragBlock on alloc of : "<< bytes <<endl;
			free( defragBlock );
			defragBlock = NULL;
			ptr = (sf_u8 *) malloc( bytes + 16 + 4 );			
			AllocDefragBlock();
		}
		if ( !ptr ) {
			Debug::debug(Debug::error,__FUNCTION__) <<"malloc failure for : "<< bytes << endl;
		}
	}
	alignedPtr = (sf_u8 *) ( ( (uintptr_t) ptr ) + 15 & ~15 );
	if ( alignedPtr - ptr < 4 ) {
		alignedPtr += 16;
	}
	*((int *)(alignedPtr - 4)) = (uintptr_t) ptr;
	return (void *) alignedPtr;
}


void CHeap::Free16( void *p ) {
	free( (void *) *((int *) (( (sf_u8 *) p ) - 4)) );
}


sf_u32 CHeap::Msize( void *p ) {

	if ( !p ) {
		return 0;
	}

#if USE_LIBC_MALLOC
	#ifdef _WIN32
		return _msize( p );
	#else
		return 0;
	#endif
#else
	switch( ((sf_u8 *)(p))[-1] ) {
		case SMALL_ALLOC: {
			return SMALL_ALIGN( ((sf_u8 *)(p))[-SMALL_HEADER_SIZE] * HEAP_ALIGN );
		}
		case MEDIUM_ALLOC: {
			return ((mediumHeapEntry_s *)(((sf_u8 *)(p)) - ALIGN_SIZE( MEDIUM_HEADER_SIZE )))->size - ALIGN_SIZE( MEDIUM_HEADER_SIZE );
		}
		case LARGE_ALLOC: {
			return ((CHeap::page_s*)(*((sf_u32 *)(((sf_u8 *)p) - ALIGN_SIZE( LARGE_HEADER_SIZE )))))->dataSize - ALIGN_SIZE( LARGE_HEADER_SIZE );
		}
		default: {
			//Todo: Debug::debug(Debug::error,__FUNCTION__) << "CHeap::Msize: invalid memory block -> "<< GetCallStackCurStr( 4 ) << endl;
			return 0;
		}
	}
#endif
}


void CHeap::Dump() {
	CHeap::page_s	*pg;

	for ( pg = smallFirstUsedPage; pg; pg = pg->next ) {
		Debug::debug(Debug::heap,__FUNCTION__) << pg->data << "  bytes "<< pg->dataSize <<"-8d  (in use by small heap)"<<endl;
	}

	if ( smallCurPage ) {
		pg = smallCurPage;
		Debug::debug(Debug::heap,__FUNCTION__) << pg->data << "  bytes "<<pg->dataSize <<"-8d  (small heap active page)"<< endl;
	}

	for ( pg = mediumFirstUsedPage; pg; pg = pg->next ) {
		Debug::debug(Debug::heap,__FUNCTION__) << pg->data<<"  bytes "<< pg->dataSize<<"-8d  (completely used by medium heap)\n"<< endl;
	}

	for ( pg = mediumFirstFreePage; pg; pg = pg->next ) {
		Debug::debug(Debug::heap,__FUNCTION__) << pg->data<<"  bytes "<< pg->dataSize <<"-8d  (partially used by medium heap)\n"<< endl;
	}
	
	for ( pg = largeFirstUsedPage; pg; pg = pg->next ) {
		Debug::debug(Debug::heap,__FUNCTION__) << pg->data<<"  bytes "<<pg->dataSize <<"-8d  (fully used by large heap)\n"<< endl;
	}

	Debug::debug(Debug::heap,__FUNCTION__) << "pages allocated : "<< pagesAllocated <<endl;
}


void CHeap::FreePageReal( CHeap::page_s *p ) {
	SMF_ASSERT( p );
	::free( p );
}


void CHeap::ReleaseSwappedPages () {
	if ( swapPage ) {
		FreePageReal( swapPage );
	}
	swapPage = NULL;
}


CHeap::page_s* CHeap::AllocatePage( sf_u32 bytes ) {
	CHeap::page_s*	p;

	pageRequests++;

	if ( swapPage && swapPage->dataSize == bytes ) {			// if we've got a swap page somewhere
		p			= swapPage;
		swapPage	= NULL;
	}
	else {
		sf_u32 size;

		size = bytes + sizeof(CHeap::page_s);

		p = (CHeap::page_s *) ::malloc( size + HEAP_ALIGN - 1 );
		if ( !p ) {
			if ( defragBlock ) {
				Debug::debug(Debug::heap,__FUNCTION__) <<  "Freeing defragBlock on alloc of : "<< size + HEAP_ALIGN - 1 <<endl;
				free( defragBlock );
				defragBlock = NULL;
				p = (CHeap::page_s *) ::malloc( size + HEAP_ALIGN - 1 );			
				AllocDefragBlock();
			}
			if ( !p ) {
				Debug::debug(Debug::error,__FUNCTION__) <<  "malloc failure for "<< bytes <<endl;
			}
		}

		p->data		= (void *) ALIGN_SIZE( (uintptr_t)((sf_u8 *)(p)) + sizeof( CHeap::page_s ) );
		p->dataSize	= size - sizeof(CHeap::page_s);
		p->firstFree = NULL;
		p->largestFree = 0;
		OSAllocs++;
	}

	p->prev = NULL;
	p->next = NULL;

	pagesAllocated++;
	
	return p;
}


void CHeap::FreePage( CHeap::page_s *p ) {
	SMF_ASSERT( p );

	if ( p->dataSize == pageSize && !swapPage ) {			// add to swap list?
		swapPage = p;
	}
	else {
		FreePageReal( p );
	}

	pagesAllocated--;
}


void *CHeap::SmallAllocate( sf_u32 bytes ) {
	// we need the at least sizeof( sf_u32 ) bytes for the free list
	if ( bytes < sizeof( sf_u32 ) ) {
		bytes = sizeof( sf_u32 );
	}

	// increase the number of bytes if necessary to make sure the next small allocation is aligned
	bytes = SMALL_ALIGN( bytes );

	sf_u8 *smallBlock = (sf_u8 *)(smallFirstFree[bytes / HEAP_ALIGN]);
	if ( smallBlock ) {
		sf_u32 *link = (sf_u32 *)(smallBlock + SMALL_HEADER_SIZE);
		smallBlock[1] = SMALL_ALLOC;					// allocation identifier
		smallFirstFree[bytes / HEAP_ALIGN] = (void *)(*link);
		return (void *)(link);
	}

	sf_u32 bytesLeft = (long)(pageSize) - smallCurPageOffset;
	// if we need to allocate a new page
	if ( bytes >= bytesLeft ) {

		smallCurPage->next	= smallFirstUsedPage;
		smallFirstUsedPage	= smallCurPage;
		smallCurPage		= AllocatePage( pageSize );
		if ( !smallCurPage ) {
			return NULL;
		}
		// make sure the first allocation is aligned
		smallCurPageOffset	= SMALL_ALIGN( 0 );
	}

	smallBlock			= ((sf_u8 *)smallCurPage->data) + smallCurPageOffset;
	smallBlock[0]		= (sf_u8)(bytes / HEAP_ALIGN);		// write # of bytes/ALIGN
	smallBlock[1]		= SMALL_ALLOC;					// allocation identifier
	smallCurPageOffset  += bytes + SMALL_HEADER_SIZE;	// increase the offset on the current page
	return ( smallBlock + SMALL_HEADER_SIZE );			// skip the first two bytes
}


void CHeap::SmallFree( void *ptr ) {
	((sf_u8 *)(ptr))[-1] = INVALID_ALLOC;

	sf_u8 *d = ( (sf_u8 *)ptr ) - SMALL_HEADER_SIZE;
	sf_u32 *dt = (sf_u32 *)ptr;
	// index into the table with free small memory blocks
	sf_u32 ix = *d;

	// check if the index is correct
	if ( ix > (256 / HEAP_ALIGN) ) {
		Debug::debug(Debug::error,__FUNCTION__) << "SmallFree: invalid memory block" <<endl;
	}

	*dt = (uintptr_t)smallFirstFree[ix];	// write next index
	smallFirstFree[ix] = (void *)d;		// link
}


void *CHeap::MediumAllocateFromPage( CHeap::page_s *p, sf_u32 sizeNeeded ) {

	mediumHeapEntry_s	*best,*nw = NULL;
	sf_u8				*ret;

	best = (mediumHeapEntry_s *)(p->firstFree);			// first block is largest

	SMF_ASSERT( best );
	SMF_ASSERT( best->size == p->largestFree );
	SMF_ASSERT( best->size >= sizeNeeded );

	// if we can allocate another block from this page after allocating sizeNeeded bytes
	if ( best->size >= (sf_u32)( sizeNeeded + MEDIUM_SMALLEST_SIZE ) ) {
		nw = (mediumHeapEntry_s *)((sf_u8 *)best + best->size - sizeNeeded);
		nw->page		= p;
		nw->prev		= best;
		nw->next		= best->next;
		nw->prevFree	= NULL;
		nw->nextFree	= NULL;
		nw->size		= sizeNeeded;
		nw->freeBlock	= 0;			// used block
		if ( best->next ) {
			best->next->prev = nw;
		}
		best->next	= nw;
		best->size	-= sizeNeeded;
		
		p->largestFree = best->size;
	}
	else {
		if ( best->prevFree ) {
			best->prevFree->nextFree = best->nextFree;
		}
		else {
			p->firstFree = (void *)best->nextFree;
		}
		if ( best->nextFree ) {
			best->nextFree->prevFree = best->prevFree;
		}

		best->prevFree  = NULL;
		best->nextFree  = NULL;
		best->freeBlock = 0;			// used block
		nw = best;

		p->largestFree = 0;
	}

	ret		= (sf_u8 *)(nw) + ALIGN_SIZE( MEDIUM_HEADER_SIZE );
	ret[-1] = MEDIUM_ALLOC;		// allocation identifier

	return (void *)(ret);
}


void *CHeap::MediumAllocate( sf_u32 bytes ) {
	CHeap::page_s		*p;
	void				*data;

	sf_u32 sizeNeeded = ALIGN_SIZE( bytes ) + ALIGN_SIZE( MEDIUM_HEADER_SIZE );

	// find first page with enough space
	for ( p = mediumFirstFreePage; p; p = p->next ) {
		if ( p->largestFree >= sizeNeeded ) {
			break;
		}
	}

	if ( !p ) {								// need to allocate new page?
		p = AllocatePage( pageSize );
		if ( !p ) {
			return NULL;					// malloc failure!
		}
		p->prev		= NULL;
		p->next		= mediumFirstFreePage;
		if (p->next) {
			p->next->prev = p;
		}
		else {
			mediumLastFreePage	= p;
		}

		mediumFirstFreePage		= p;
		
		p->largestFree	= pageSize;
		p->firstFree	= (void *)p->data;

		mediumHeapEntry_s *e;
		e				= (mediumHeapEntry_s *)(p->firstFree);
		e->page			= p;
		// make sure ((sf_u8 *)e + e->size) is aligned
		e->size			= pageSize & ~(HEAP_ALIGN - 1);
		e->prev			= NULL;
		e->next			= NULL;
		e->prevFree		= NULL;
		e->nextFree		= NULL;
		e->freeBlock	= 1;
	}

	data = MediumAllocateFromPage( p, sizeNeeded );		// allocate data from page

    // if the page can no longer serve memory, move it away from free list
	// (so that it won't slow down the later alloc queries)
	// this modification speeds up the pageWalk from O(N) to O(sqrt(N))
	// a call to free may swap this page back to the free list

	if ( p->largestFree < MEDIUM_SMALLEST_SIZE ) {
		if ( p == mediumLastFreePage ) {
			mediumLastFreePage = p->prev;
		}

		if ( p == mediumFirstFreePage ) {
			mediumFirstFreePage = p->next;
		}

		if ( p->prev ) {
			p->prev->next = p->next;
		}
		if ( p->next ) {
			p->next->prev = p->prev;
		}

		// link to "completely used" list
		p->prev = NULL;
		p->next = mediumFirstUsedPage;
		if ( p->next ) {
			p->next->prev = p;
		}
		mediumFirstUsedPage = p;
		return data;
	} 

	// re-order linked list (so that next malloc query starts from current
	// matching block) -- this speeds up both the page walks and block walks

	if ( p != mediumFirstFreePage ) {
		SMF_ASSERT( mediumLastFreePage );
		SMF_ASSERT( mediumFirstFreePage );
		SMF_ASSERT( p->prev);

		mediumLastFreePage->next	= mediumFirstFreePage;
		mediumFirstFreePage->prev	= mediumLastFreePage;
		mediumLastFreePage			= p->prev;
		p->prev->next				= NULL;
		p->prev						= NULL;
		mediumFirstFreePage			= p;
	}

	return data;
}


void CHeap::MediumFree( void *ptr ) {
	((sf_u8 *)(ptr))[-1] = INVALID_ALLOC;

	mediumHeapEntry_s	*e = (mediumHeapEntry_s *)((sf_u8 *)ptr - ALIGN_SIZE( MEDIUM_HEADER_SIZE ));
	CHeap::page_s		*p = e->page;
	bool				isInFreeList;

	isInFreeList = p->largestFree >= MEDIUM_SMALLEST_SIZE;

	SMF_ASSERT( e->size );
	SMF_ASSERT( e->freeBlock == 0 );

	mediumHeapEntry_s *prev = e->prev;

	// if the previous block is free we can merge
	if ( prev && prev->freeBlock ) {
		prev->size += e->size;
		prev->next = e->next;
		if ( e->next ) {
			e->next->prev = prev;
		}
		e = prev;
	}
	else {
		e->prevFree		= NULL;				// link to beginning of free list
		e->nextFree		= (mediumHeapEntry_s *)p->firstFree;
		if ( e->nextFree ) {
			SMF_ASSERT( !(e->nextFree->prevFree) );
			e->nextFree->prevFree = e;
		}

		p->firstFree	= e;
		p->largestFree	= e->size;
		e->freeBlock	= 1;				// mark block as free
	}
			
	mediumHeapEntry_s *next = e->next;

	// if the next block is free we can merge
	if ( next && next->freeBlock ) {
		e->size += next->size;
		e->next = next->next;
		
		if ( next->next ) {
			next->next->prev = e;
		}
		
		if ( next->prevFree ) {
			next->prevFree->nextFree = next->nextFree;
		}
		else {
			SMF_ASSERT( next == p->firstFree );
			p->firstFree = next->nextFree;
		}

		if ( next->nextFree ) {
			next->nextFree->prevFree = next->prevFree;
		}
	}

	if ( p->firstFree ) {
		p->largestFree = ((mediumHeapEntry_s *)(p->firstFree))->size;
	}
	else {
		p->largestFree = 0;
	}

	// did e become the largest block of the page ?

	if ( e->size > p->largestFree ) {
		SMF_ASSERT( e != p->firstFree );
		p->largestFree = e->size;

		if ( e->prevFree ) {
			e->prevFree->nextFree = e->nextFree;
		}
		if ( e->nextFree ) {
			e->nextFree->prevFree = e->prevFree;
		}
		
		e->nextFree = (mediumHeapEntry_s *)p->firstFree;
		e->prevFree = NULL;
		if ( e->nextFree ) {
			e->nextFree->prevFree = e;
		}
		p->firstFree = e;
	}

	// if page wasn't in free list (because it was near-full), move it back there
	if ( !isInFreeList ) {

		// remove from "completely used" list
		if ( p->prev ) {
			p->prev->next = p->next;
		}
		if ( p->next ) {
			p->next->prev = p->prev;
		}
		if ( p == mediumFirstUsedPage ) {
			mediumFirstUsedPage = p->next;
		}

		p->next = NULL;
		p->prev = mediumLastFreePage;

		if ( mediumLastFreePage ) {
			mediumLastFreePage->next = p;
		}
		mediumLastFreePage = p;
		if ( !mediumFirstFreePage ) {
			mediumFirstFreePage = p;
		}
	} 
}


void *CHeap::LargeAllocate( sf_u32 bytes ) {
	CHeap::page_s *p = AllocatePage( bytes + ALIGN_SIZE( LARGE_HEADER_SIZE ) );

	SMF_ASSERT( p );

	if ( !p ) {
		return NULL;
	}

	sf_u8 *	d	= (sf_u8*)(p->data) + ALIGN_SIZE( LARGE_HEADER_SIZE );
	sf_u32 *	dw	= (sf_u32*)(d - ALIGN_SIZE( LARGE_HEADER_SIZE ));
	dw[0]		= (uintptr_t)p;				// write pointer back to page table
	d[-1]		= LARGE_ALLOC;			// allocation identifier

	// link to 'large used page list'
	p->prev = NULL;
	p->next = largeFirstUsedPage;
	if ( p->next ) {
		p->next->prev = p;
	}
	largeFirstUsedPage = p;

	return (void *)(d);
}


void CHeap::LargeFree( void *ptr) {
	CHeap::page_s*	pg;

	((sf_u8 *)(ptr))[-1] = INVALID_ALLOC;

	// get page pointer
	pg = (CHeap::page_s *)(*((sf_u32 *)(((sf_u8 *)ptr) - ALIGN_SIZE( LARGE_HEADER_SIZE ))));

	// unlink from doubly linked list
	if ( pg->prev ) {
		pg->prev->next = pg->next;
	}
	if ( pg->next ) {
		pg->next->prev = pg->prev;
	}
	if ( pg == largeFirstUsedPage ) {
		largeFirstUsedPage = pg->next;
	}
	pg->next = pg->prev = NULL;

	FreePage(pg);
}

//===============================================================
//
//	memory allocation all in one place
//
//===============================================================

#undef new

static CHeap *			mem_heap = NULL;
static memoryStats_t	mem_total_allocs = { 0, 0x0fffffff, -1, 0 };
static memoryStats_t	mem_frame_allocs;
static memoryStats_t	mem_frame_frees;


void mem_ClearFrameStats() {
	mem_frame_allocs.num = mem_frame_frees.num = 0;
	mem_frame_allocs.minSize = mem_frame_frees.minSize = 0x0fffffff;
	mem_frame_allocs.maxSize = mem_frame_frees.maxSize = -1;
	mem_frame_allocs.totalSize = mem_frame_frees.totalSize = 0;
}


void mem_GetFrameStats( memoryStats_t &allocs, memoryStats_t &frees ) {
	allocs = mem_frame_allocs;
	frees = mem_frame_frees;
}


void mem_GetStats( memoryStats_t &stats ) {
	stats = mem_total_allocs;
}


void mem_UpdateStats( memoryStats_t &stats, int size ) {
	stats.num++;
	if ( size < stats.minSize ) {
		stats.minSize = size;
	}
	if ( size > stats.maxSize ) {
		stats.maxSize = size;
	}
	stats.totalSize += size;
}


void mem_UpdateAllocStats( int size ) {
	mem_UpdateStats( mem_frame_allocs, size );
	mem_UpdateStats( mem_total_allocs, size );
}


void mem_UpdateFreeStats( int size ) {
	mem_UpdateStats( mem_frame_frees, size );
	mem_total_allocs.num--;
	mem_total_allocs.totalSize -= size;
}


#ifndef SMF_DEBUG_MEM

void *mem_Alloc( const int size ) {
	if ( !size ) {
		return NULL;
	}
	if ( !mem_heap ) {
#ifdef CRASH_ON_STATIC_ALLOCATION
		*((int*)0x0) = 1;
#endif
		return malloc( size );
	}
	void *mem = mem_heap->Allocate( size );
	mem_UpdateAllocStats( mem_heap->Msize( mem ) );
	return mem;
}


void mem_Free( void *ptr ) {
	if ( !ptr ) {
		return;
	}
	if ( !mem_heap ) {
#ifdef CRASH_ON_STATIC_ALLOCATION
		*((int*)0x0) = 1;
#endif
		free( ptr );
		return;
	}
	mem_UpdateFreeStats( mem_heap->Msize( ptr ) );
 	mem_heap->free( ptr );
}


void *mem_Alloc16( const int size ) {
	if ( !size ) {
		return NULL;
	}
	if ( !mem_heap ) {
#ifdef CRASH_ON_STATIC_ALLOCATION
		*((int*)0x0) = 1;
#endif
		return malloc( size );
	}
	void *mem = mem_heap->Allocate16( size );
	// make sure the memory is 16 byte aligned
	SMF_ASSERT( ( ((uintptr_t)mem) & 15) == 0 );
	return mem;
}


void mem_Free16( void *ptr ) {
	if ( !ptr ) {
		return;
	}
	if ( !mem_heap ) {
#ifdef CRASH_ON_STATIC_ALLOCATION
		*((int*)0x0) = 1;
#endif
		free( ptr );
		return;
	}
	// make sure the memory is 16 sf_u8 aligned
	SMF_ASSERT( ( ((uintptr_t)ptr) & 15) == 0 );
 	mem_heap->Free16( ptr );
}


void *mem_ClearedAlloc( const int size ) {
	void *mem = mem_Alloc( size );
	MATH::SIMDProcessor->memSet( mem, 0, size );
	return mem;
}

/*
==================
mem_ClearedAlloc
==================
*/
void mem_AllocDefragBlock() {
	mem_heap->AllocDefragBlock();
}


char *mem_CopyString( const char *in ) {
	char	*out;
	
	out = (char *)mem_Alloc( strlen(in) + 1 );
	strcpy( out, in );
	return out;
}


void mem_Dump_f( const CCMDLineArgs &args ) {
}


void mem_DumpCompressed_f( const CCMDLineArgs &args ) {
}


void mem_Init() {
	mem_heap = new CHeap;
	mem_ClearFrameStats();
}
bool mem_isInitialized(){
	return mem_heap->isInitialized();
}


void mem_Shutdown() {
	CHeap *m = mem_heap;
	mem_heap = NULL;
	delete m;
}


void mem_EnableLeakTest( const char *name ) {
}


#else /* !SMF_DEBUG_MEM */

#undef		mem_Alloc
#undef		mem_ClearedAlloc
#undef		Com_ClearedReAlloc
#undef		mem_Free
#undef		mem_CopyString
#undef		mem_Alloc16
#undef		mem_Free16

#define MAX_CALLSTACK_DEPTH		6

// size of this struct must be a multiple of 16 bytes
typedef struct debugMemory_s {
	const char *			fileName;
	int						lineNumber;
	int						frameNumber;
	int						size;
	address_t				callStack[MAX_CALLSTACK_DEPTH];
	struct debugMemory_s *	prev;
	struct debugMemory_s *	next;
} debugMemory_t;

static debugMemory_t *	mem_debugMemory = NULL;
static char				mem_leakName[256] = "";

/*
==================
mem_CleanupFileName
==================
*/
const char *mem_CleanupFileName( const char *fileName ) {
	int i1, i2;
	CMyString newFileName;
	static char newFileNames[4][MAX_STRING_CHARS];
	static int index;

	newFileName = fileName;
	newFileName.backSlashesToSlashes();
	i1 = newFileName.find( "neo", false );
	if ( i1 >= 0 ) {
		i1 = newFileName.find( "/", false, i1 );
		newFileName = newFileName.right( newFileName.getLenght() - ( i1 + 1 ) );
	}
	while( 1 ) {
		i1 = newFileName.find( "/../" );
		if ( i1 <= 0 ) {
			break;
		}
		i2 = i1 - 1;
		while( i2 > 1 && newFileName[i2-1] != '/' ) {
			i2--;
		}
		newFileName = newFileName.left( i2 - 1 ) + newFileName.right( newFileName.getLenght() - ( i1 + 4 ) );
	}
	index = ( index + 1 ) & 3;
	strncpy( newFileNames[index], newFileName.c_str(), sizeof( newFileNames[index] ) );
	return newFileNames[index];
}

/*
==================
mem_Dump
==================
*/
void mem_Dump( const char *fileName ) {
	int i, numBlocks, totalSize;
	char dump[32], *ptr;
	debugMemory_t *b;
	CMyString module, funcName;
	FILE *f;

	f = fopen( fileName, "wb" );
	if ( !f ) {
		return;
	}

	totalSize = 0;
	for ( numBlocks = 0, b = mem_debugMemory; b; b = b->next, numBlocks++ ) {
		ptr = ((char *) b) + sizeof(debugMemory_t);
		totalSize += b->size;
		for ( i = 0; i < (sizeof(dump)-1) && i < b->size; i++) {
			if ( ptr[i] >= 32 && ptr[i] < 127 ) {
				dump[i] = ptr[i];
			} else {
				dump[i] = '_';
			}
		}
		dump[i] = '\0';
		if ( ( b->size >> 10 ) != 0 ) {
			fprintf( f, "size: %6d KB: %s, line: %d [%s], call stack: %s\r\n", ( b->size >> 10 ), mem_CleanupFileName(b->fileName), b->lineNumber, dump, idLib::sys->GetCallStackStr( b->callStack, MAX_CALLSTACK_DEPTH ) );
		}
		else {
			fprintf( f, "size: %7d B: %s, line: %d [%s], call stack: %s\r\n", b->size, mem_CleanupFileName(b->fileName), b->lineNumber, dump, idLib::sys->GetCallStackStr( b->callStack, MAX_CALLSTACK_DEPTH ) );
		}
	}

	idLib::sys->ShutdownSymbols();

	fprintf( f, "%8d total memory blocks allocated\r\n", numBlocks );
	fprintf( f, "%8d KB memory allocated\r\n", ( totalSize >> 10 ) );

	fclose( f );
}

/*
==================
mem_Dump_f
==================
*/
void mem_Dump_f( const CCMDLineArgs &args ) {
	const char *fileName;

	if ( args.Argc() >= 2 ) {
		fileName = args.Argv( 1 );
	}
	else {
		fileName = "memorydump.txt";
	}
	mem_Dump( fileName );
}

/*
==================
mem_DumpCompressed
==================
*/
typedef struct allocInfo_s {
	const char *			fileName;
	int						lineNumber;
	int						size;
	int						numAllocs;
	address_t				callStack[MAX_CALLSTACK_DEPTH];
	struct allocInfo_s *	next;
} allocInfo_t;

typedef enum {
	MEMSORT_SIZE,
	MEMSORT_LOCATION,
	MEMSORT_NUMALLOCS,
	MEMSORT_CALLSTACK
} memorySortType_t;

void mem_DumpCompressed( const char *fileName, memorySortType_t memSort, int sortCallStack, int numFrames ) {
	int numBlocks, totalSize, r, j;
	debugMemory_t *b;
	allocInfo_t *a, *nexta, *allocInfo = NULL, *sortedAllocInfo = NULL, *prevSorted, *nextSorted;
	CMyString module, funcName;
	FILE *f;

	// build list with memory allocations
	totalSize = 0;
	numBlocks = 0;
	for ( b = mem_debugMemory; b; b = b->next ) {

		if ( numFrames && b->frameNumber < idLib::frameNumber - numFrames ) {
			continue;
		}

		numBlocks++;
		totalSize += b->size;

		// search for an allocation from the same source location
		for ( a = allocInfo; a; a = a->next ) {
			if ( a->lineNumber != b->lineNumber ) {
				continue;
			}
			for ( j = 0; j < MAX_CALLSTACK_DEPTH; j++ ) {
				if ( a->callStack[j] != b->callStack[j] ) {
					break;
				}
			}
			if ( j < MAX_CALLSTACK_DEPTH ) {
				continue;
			}
			if ( CMyString::compare( a->fileName, b->fileName ) != 0 ) {
				continue;
			}
			a->numAllocs++;
			a->size += b->size;
			break;
		}

		// if this is an allocation from a new source location
		if ( !a ) {
			a = (allocInfo_t *) ::malloc( sizeof( allocInfo_t ) );
			a->fileName = b->fileName;
			a->lineNumber = b->lineNumber;
			a->size = b->size;
			a->numAllocs = 1;
			for ( j = 0; j < MAX_CALLSTACK_DEPTH; j++ ) {
				a->callStack[j] = b->callStack[j];
			}
			a->next = allocInfo;
			allocInfo = a;
		}
	}

	// sort list
	for ( a = allocInfo; a; a = nexta ) {
		nexta = a->next;

		prevSorted = NULL;
		switch( memSort ) {
			// sort on size
			case MEMSORT_SIZE: {
				for ( nextSorted = sortedAllocInfo; nextSorted; nextSorted = nextSorted->next ) {
					if ( a->size > nextSorted->size ) {
						break;
					}
					prevSorted = nextSorted;
				}
				break;
			}
			// sort on file name and line number
			case MEMSORT_LOCATION: {
				for ( nextSorted = sortedAllocInfo; nextSorted; nextSorted = nextSorted->next ) {
					r = CMyString::compare( mem_CleanupFileName( a->fileName ), mem_CleanupFileName( nextSorted->fileName ) );
					if ( r < 0 || ( r == 0 && a->lineNumber < nextSorted->lineNumber ) ) {
						break;
					}
					prevSorted = nextSorted;
				}
				break;
			}
			// sort on the number of allocations
			case MEMSORT_NUMALLOCS: {
				for ( nextSorted = sortedAllocInfo; nextSorted; nextSorted = nextSorted->next ) {
					if ( a->numAllocs > nextSorted->numAllocs ) {
						break;
					}
					prevSorted = nextSorted;
				}
				break;
			}
			// sort on call stack
			case MEMSORT_CALLSTACK: {
				for ( nextSorted = sortedAllocInfo; nextSorted; nextSorted = nextSorted->next ) {
					if ( a->callStack[sortCallStack] < nextSorted->callStack[sortCallStack] ) {
						break;
					}
					prevSorted = nextSorted;
				}
				break;
			}
		}
		if ( !prevSorted ) {
			a->next = sortedAllocInfo;
			sortedAllocInfo = a;
		}
		else {
			prevSorted->next = a;
			a->next = nextSorted;
		}
	}

	f = fopen( fileName, "wb" );
	if ( !f ) {
		return;
	}

	// write list to file
	for ( a = sortedAllocInfo; a; a = nexta ) {
		nexta = a->next;
		fprintf( f, "size: %6d KB, allocs: %5d: %s, line: %d, call stack: %s\r\n",
					(a->size >> 10), a->numAllocs, mem_CleanupFileName(a->fileName),
							a->lineNumber, idLib::sys->GetCallStackStr( a->callStack, MAX_CALLSTACK_DEPTH ) );
		::free( a );
	}

	idLib::sys->ShutdownSymbols();

	fprintf( f, "%8d total memory blocks allocated\r\n", numBlocks );
	fprintf( f, "%8d KB memory allocated\r\n", ( totalSize >> 10 ) );

	fclose( f );
}

/*
==================
mem_DumpCompressed_f
==================
*/
void mem_DumpCompressed_f( const CCMDLineArgs &args ) {
	int argNum;
	const char *arg, *fileName;
	memorySortType_t memSort = MEMSORT_LOCATION;
	int sortCallStack = 0, numFrames = 0;

	// get cmd-line options
	argNum = 1;
	arg = args.Argv( argNum );
	while( arg[0] == '-' ) {
		arg = args.Argv( ++argNum );
		if ( CMyString::compareInsen( arg, "s" ) == 0 ) {
			memSort = MEMSORT_SIZE;
		} else if ( CMyString::compareInsen( arg, "l" ) == 0 ) {
			memSort = MEMSORT_LOCATION;
		} else if ( CMyString::compareInsen( arg, "a" ) == 0 ) {
			memSort = MEMSORT_NUMALLOCS;
		} else if ( CMyString::compareInsen( arg, "cs1" ) == 0 ) {
			memSort = MEMSORT_CALLSTACK;
			sortCallStack = 2;
		} else if ( CMyString::compareInsen( arg, "cs2" ) == 0 ) {
			memSort = MEMSORT_CALLSTACK;
			sortCallStack = 1;
		} else if ( CMyString::compareInsen( arg, "cs3" ) == 0 ) {
			memSort = MEMSORT_CALLSTACK;
			sortCallStack = 0;
		} else if ( arg[0] == 'f' ) {
			numFrames = atoi( arg + 1 );
		} else {
			Debug::debug(Debug::math,__FUNCTION__) << "memoryDumpCompressed [options] [filename]\n"
						"options:\n"
						"  -s     sort on size\n"
						"  -l     sort on location\n"
						"  -a     sort on the number of allocations\n"
						"  -cs1   sort on first function on call stack\n"
						"  -cs2   sort on second function on call stack\n"
						"  -cs3   sort on third function on call stack\n"
						"  -f<X>  only report allocations the last X frames\n"
						"By default the memory allocations are sorted on location.\n"
						"By default a 'memorydump.txt' is written if no file name is specified.\n" );
			return;
		}
		arg = args.Argv( ++argNum );
	}
	if ( argNum >= args.Argc() ) {
		fileName = "memorydump.txt";
	} else {
		fileName = arg;
	}
	mem_DumpCompressed( fileName, memSort, sortCallStack, numFrames );
}

/*
==================
mem_AllocDebugMemory
==================
*/
void *mem_AllocDebugMemory( const int size, const char *fileName, const int lineNumber, const bool align16 ) {
	void *p;
	debugMemory_t *m;

	if ( !size ) {
		return NULL;
	}

	if ( !mem_heap ) {
#ifdef CRASH_ON_STATIC_ALLOCATION
		*((int*)0x0) = 1;
#endif
		// NOTE: set a breakpoint here to find memory allocations before mem_heap is initialized
		return malloc( size );
	}

	if ( align16 ) {
		p = mem_heap->Allocate16( size + sizeof( debugMemory_t ) );
	}
	else {
		p = mem_heap->Allocate( size + sizeof( debugMemory_t ) );
	}

	mem_UpdateAllocStats( size );

	m = (debugMemory_t *) p;
	m->fileName = fileName;
	m->lineNumber = lineNumber;
	m->frameNumber = idLib::frameNumber;
	m->size = size;
	m->next = mem_debugMemory;
	m->prev = NULL;
	if ( mem_debugMemory ) {
		mem_debugMemory->prev = m;
	}
	mem_debugMemory = m;
	idLib::sys->GetCallStack( m->callStack, MAX_CALLSTACK_DEPTH );

	return ( ( (sf_u8 *) p ) + sizeof( debugMemory_t ) );
}

/*
==================
mem_FreeDebugMemory
==================
*/
void mem_FreeDebugMemory( void *p, const char *fileName, const int lineNumber, const bool align16 ) {
	debugMemory_t *m;

	if ( !p ) {
		return;
	}

	if ( !mem_heap ) {
#ifdef CRASH_ON_STATIC_ALLOCATION
		*((int*)0x0) = 1;
#endif
		// NOTE: set a breakpoint here to find memory being freed before mem_heap is initialized
		free( p );
		return;
	}

	m = (debugMemory_t *) ( ( (sf_u8 *) p ) - sizeof( debugMemory_t ) );

	if ( m->size < 0 ) {
		idLib::common->FatalError( "memory freed twice, first from %s, now from %s", idLib::sys->GetCallStackStr( m->callStack, MAX_CALLSTACK_DEPTH ), idLib::sys->GetCallStackCurStr( MAX_CALLSTACK_DEPTH ) );
	}

	mem_UpdateFreeStats( m->size );

	if ( m->next ) {
		m->next->prev = m->prev;
	}
	if ( m->prev ) {
		m->prev->next = m->next;
	}
	else {
		mem_debugMemory = m->next;
	}

	m->fileName = fileName;
	m->lineNumber = lineNumber;
	m->frameNumber = idLib::frameNumber;
	m->size = -m->size;
	idLib::sys->GetCallStack( m->callStack, MAX_CALLSTACK_DEPTH );

	if ( align16 ) {
 		mem_heap->Free16( m );
	}
	else {
 		mem_heap->free( m );
	}
}

/*
==================
mem_Alloc
==================
*/
void *mem_Alloc( const int size, const char *fileName, const int lineNumber ) {
	if ( !size ) {
		return NULL;
	}
	return mem_AllocDebugMemory( size, fileName, lineNumber, false );
}

/*
==================
mem_Free
==================
*/
void mem_Free( void *ptr, const char *fileName, const int lineNumber ) {
	if ( !ptr ) {
		return;
	}
	mem_FreeDebugMemory( ptr, fileName, lineNumber, false );
}

/*
==================
mem_Alloc16
==================
*/
void *mem_Alloc16( const int size, const char *fileName, const int lineNumber ) {
	if ( !size ) {
		return NULL;
	}
	void *mem = mem_AllocDebugMemory( size, fileName, lineNumber, true );
	// make sure the memory is 16 sf_u8 aligned
	SMF_ASSERT( ( ((int)mem) & 15) == 0 );
	return mem;
}

/*
==================
mem_Free16
==================
*/
void mem_Free16( void *ptr, const char *fileName, const int lineNumber ) {
	if ( !ptr ) {
		return;
	}
	// make sure the memory is 16 sf_u8 aligned
	SMF_ASSERT( ( ((int)ptr) & 15) == 0 );
	mem_FreeDebugMemory( ptr, fileName, lineNumber, true );
}

/*
==================
mem_ClearedAlloc
==================
*/
void *mem_ClearedAlloc( const int size, const char *fileName, const int lineNumber ) {
	void *mem = mem_Alloc( size, fileName, lineNumber );
	SIMDProcessor->memSet( mem, 0, size );
	return mem;
}

/*
==================
mem_CopyString
==================
*/
char *mem_CopyString( const char *in, const char *fileName, const int lineNumber ) {
	char	*out;
	
	out = (char *)mem_Alloc( strlen(in) + 1, fileName, lineNumber );
	strcpy( out, in );
	return out;
}

/*
==================
mem_Init
==================
*/
void mem_Init() {
	mem_heap = new CHeap;
}

/*
==================
mem_Shutdown
==================
*/
void mem_Shutdown() {

	if ( mem_leakName[0] != '\0' ) {
		mem_DumpCompressed( va( "%s_leak_size.txt", mem_leakName ), MEMSORT_SIZE, 0, 0 );
		mem_DumpCompressed( va( "%s_leak_location.txt", mem_leakName ), MEMSORT_LOCATION, 0, 0 );
		mem_DumpCompressed( va( "%s_leak_cs1.txt", mem_leakName ), MEMSORT_CALLSTACK, 2, 0 );
	}

	CHeap *m = mem_heap;
	mem_heap = NULL;
	delete m;
}

/*
==================
mem_EnableLeakTest
==================
*/
void mem_EnableLeakTest( const char *name ) {
	CMyString::copyNz( mem_leakName, name, sizeof( mem_leakName ) );
}

#endif /* !SMF_DEBUG_MEM */


void * mem_alignedMalloc(size_t size, size_t alignment)
{
	alignment = alignment > 0 ? alignment : 1;
	size += alignment;

	uintptr_t ptr = (uintptr_t)malloc(size);
#ifndef __clang_analyzer__ // Hide clang analyzer false positive: Memory is never released; potential leak of memory pointed to by 'ptr'
	if (!ptr)
		return 0;
	++ptr; // Must make room for storing the offset info.
	ptrdiff_t incr = (alignment - (ptr & (alignment-1))) & (alignment-1);
	ptr += incr;
	((sf_u8*)ptr)[-1] = (sf_u8)(incr+1);
#endif
	//check alignment
	assert_X_byte_aligned(ptr,alignment);
	//SMF_ASSERT(ptr % alignment == 0);
	return (void*)ptr;
}

void mem_alignedFree(void *ptr)
{
	if (!ptr)
		return;
	sf_u8 *p = (sf_u8*)ptr;
	p -= p[-1];
	free(p);
}



} //end Util
} //end SMF
