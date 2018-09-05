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


#ifndef _SMF__LIST_H__
#define _SMF__LIST_H__

#include "../SMF_Config.h"
namespace SMF{
/**
 * \class CList
 *
 * \ingroup SMF_Data_Structures
 *
 * \brief  
 * \if pt_br
 * \brief Implementa uma Lista
 * \note não é uma lista ligada. É uma lista vetorizada
 * \note ItemType é uma estrutura ou classe, para a qual se deseja criar a fila
 * \elseif us_en
 * \brief Implements a List
 * \endif
 *
 * 
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */



/*

CListSortCompare<ItemType>

*/
#ifdef __INTEL_COMPILER
// intel compiler messes thing here
template< class ItemType >
SMF_INLINE_FORCED int CListSortCompare( const ItemType *a, const ItemType *b ) {
	SMF_ASSERT( 0 );
	return 0;
}
#else
template< class ItemType >
SMF_INLINE_FORCED int CListSortCompare( const ItemType *a, const ItemType *b ) {
	return *a - *b;
}
#endif

/*

CListNewElement<ItemType>

*/

template< class ItemType >
SMF_INLINE_FORCED ItemType *CListNewElement() {
	return new ItemType;
}

/*

swapElements<ItemType>

*/
template< class ItemType >
SMF_INLINE_FORCED void swapElements( ItemType &a, ItemType &b ) {
	ItemType c = a;
	a = b;
	b = c;
}

template< class ItemType >
class CList {
public:

	typedef int		cmp_t( const ItemType *, const ItemType * );
	typedef ItemType	new_t();

					CList( int newgranularity = 16 );
					CList( const CList<ItemType> &other );
					~CList<ItemType>();
	//! limpa a lista
	// clear the list
	void			clear();				

	//! retorna o número de elementos da lista
	// returns number of elements in list
	int				getNum() const;									
	
	//! retorna o número de elementos da lista
	// returns number of elements in list
	int				getNumElements() const;	
	//! retorna o número de elementos alocadaos	
	// returns number of elements allocated for
	int				getNumAllocated() const;							
	/// set new granularity
	void			setGranularity( int newgranularity );	
	/** 
	\brief get the current granularity
	**/
	int				getGranularity() const;						
	// returns total size of allocated memory
	size_t			getAllocatedMem() const;			
	// returns total size of allocated memory including size of list ItemType
	size_t			getSize() const;			
	// returns size of the used elements in the list
	size_t			getMemoryUsed() const;							

	CList<ItemType> &	operator=( const CList<ItemType> &other );
	const ItemType &	operator[]( int index ) const;
	ItemType &			operator[]( int index );
	/// resizes list to exactly the number of elements it contains
	void			condense();
	/// resizes list to the given number of elements
	void			resize( int newsize );								
	/// resizes list and sets new granularity
	void			resize( int newsize, int newgranularity	 );			
	/// set number of elements in list and resize to exactly this number if necessary
	void			setNum( int newnum, bool resize = true );			
	/// assure list has given number of elements, but leave them uninitialized
	void			assureSize( int newSize);							
	/// assure list has given number of elements and initialize any new elements
	void			assureSize( int newSize, const ItemType &initValue );	
	/// assure the pointer list has the given number of elements and allocate any new elements
	void			assureSizeAlloc( int newSize, new_t *allocator );	
	/// returns a pointer to the list
	ItemType *			getPtr();
	/// returns a pointer to the list
	const ItemType *	getPtr() const;									
	/// returns reference to a new data element at the end of the list
	ItemType &			alloc();										
	/// append element
	int				append( const ItemType & obj );							
	/// append list
	int				append( const CList<ItemType> &other );				
	/// add unique element
	int				addUnique( const ItemType & obj );						
	/// insert the element at the given index
	int				insert( const ItemType & obj, int index = 0 );			
	/// find the index for the given element
	int				findIndex( const ItemType & obj ) const;				
	/// find pointer to the given element
	ItemType *		find( ItemType const & obj ) const;						
	/// find the index for the first NULL pointer in the list
	int				findNull() const;								
	/// returns the index for the pointer to an element in the list
	int				getIndexOf( const ItemType *obj ) const;					
	/// remove the element at the given index
	bool			removeIndex( int index );							
	/// remove the element
	bool			remove( const ItemType & obj );							
	void			sort( cmp_t *compare = ( cmp_t * )&CListSortCompare<ItemType> );
	void			sortSubSection( int startIndex, int endIndex, cmp_t *compare = ( cmp_t * )&CListSortCompare<ItemType> );
	/// swap the contents of the lists
	void			swapList( CList<ItemType> &other );						
	/// delete the contents of the list
	void			deleteContents( bool clear );						

private:
	int				m_iNum;
	int				m_iSize;
	int				m_iGranularity;
	ItemType *		m_pData;
};

/*
CList<ItemType>::CList( int )
*/
/** \brief cria uma nova lista 
* \param newgranularity: granularidade é o número de elementos cuja capacidade da lista cresce
*
*/
template< class ItemType >
SMF_INLINE_FORCED CList<ItemType>::CList( int newgranularity ) {
	SMF_ASSERT( newgranularity > 0 );

	m_pData		= NULL;
	m_iGranularity	= newgranularity;
	clear();
}

/*
CList<ItemType>::CList( const CList<ItemType> &other )
*/
/** \brief
*
*/
template< class ItemType >
SMF_INLINE_FORCED CList<ItemType>::CList( const CList<ItemType> &other ) {
	m_pData = NULL;
	*this = other;
}

/*
CList<ItemType>::~CList<ItemType>
*/
/** \brief
*
*/
template< class ItemType >
SMF_INLINE_FORCED CList<ItemType>::~CList() {
	clear();
}

/*
CList<ItemType>::clear
Frees up the memory allocated by the list.  Assumes that ItemType automatically handles freeing up memory.
*/
/** \brief
*
*/
template< class ItemType >
SMF_INLINE_FORCED void CList<ItemType>::clear() {
	if ( m_pData ) {
		delete[] m_pData;
	}

	m_pData	= NULL;
	m_iNum		= 0;
	m_iSize	= 0;
}

/*

CList<ItemType>::deleteContents

Calls the destructor of all elements in the list.  Conditionally frees up 
memory used by the list.
Note that this only works on lists containing pointers to objects and will 
cause a compiler error if called with non-pointers.  Since the list was not 
responsible for allocating the object, it has no information on whether the 
object still exists or not, so care must be taken to ensure that
the pointers are still valid when this function is called.  
Function will set all pointers in the list to NULL.

*/
/** \brief apaga o conteúdo da lista
* chama o destrutor de todos os elementos na lista. Condicionalmente  libera
* a memória consumida pela lista. 
* \note este método só funciona para lista que contém ponteiros para objetos e
* irá causar erro de compilação se utilizado com objetos não ponteiros.
* \note Deve-se tomar cuidado para garantir que os ponteiros dos objetos apontam 
* para endereços válidos quando este método for invocado, visto que a lista não 
* tem informações se os objetos ainda existem ou não, pois ela não é responsável 
* por alocá-los
* \note todos os ponteiros serão setados em NULL
*/
template< class ItemType >
SMF_INLINE_FORCED void CList<ItemType>::deleteContents( bool p_clear ) {
	int i;

	for( i = 0; i < m_iNum; i++ ) {
		delete m_pData[ i ];
		m_pData[ i ] = NULL;
	}

	if ( p_clear ) {
		clear();
	} else {
		memset( m_pData, 0, m_iSize * sizeof( ItemType ) );
	}
}

/*
CList<ItemType>::getAllocatedMem
return total memory allocated for the list in bytes, but doesn't take into account additional memory allocated by ItemType
*/
/** \brief retorna a quantidade de memória alocada
* \return size_t contendo a qantidade de memória alocada
*/
template< class ItemType >
SMF_INLINE_FORCED size_t CList<ItemType>::getAllocatedMem() const {
	return m_iSize * sizeof( ItemType );
}

/*
CList<ItemType>::getSize
return total size of list in bytes, but doesn't take into account additional memory allocated by ItemType
*/
/** \brief retorna o tamanho total da lista em bytes
* \note não leva em conta a memória adicional alocada pelo ItemType
* \return size_t contendo a quantidade de memória alocada
*/
template< class ItemType >
SMF_INLINE_FORCED size_t CList<ItemType>::getSize() const {
	return sizeof( CList<ItemType> ) + getAllocatedMem();
}

/*
CList<ItemType>::getMemoryUsed
*/
/** \brief retorna memória utilizada pelos objetos
* \return size_t contendo a quantidade de memória utilizada
*/
template< class ItemType >
SMF_INLINE_FORCED size_t CList<ItemType>::getMemoryUsed() const {
	return m_iNum * sizeof( *m_pData );
}

/*
CList<ItemType>::getNum
Returns the number of elements currently contained in the list.
Note that this is NOT an indication of the memory allocated.
*/
/** \brief retorna aquantidade de elementos contidos pela lista
* \note não é uma indicação da memória utilizada
*/
template< class ItemType >
SMF_INLINE_FORCED int CList<ItemType>::getNum() const {
	return m_iNum;
}
template< class ItemType >
SMF_INLINE_FORCED int CList<ItemType>::getNumElements() const {
	return getNum();
}
/*
CList<ItemType>::getNumAllocated
Returns the number of elements currently allocated for.
*/
/** \brief retorna o número de elementos alocados na lista
*
*/
template< class ItemType >
SMF_INLINE_FORCED int CList<ItemType>::getNumAllocated() const {
	return m_iSize;
}

/*
CList<ItemType>::setNum
resize to the exact size specified irregardless of granularity
*/
/** \brief aumenta o tamanho da lista no exato tamanho passado, 
* \note independente da granularidade
*
*/
template< class ItemType >
SMF_INLINE_FORCED void CList<ItemType>::setNum( int newnum, bool resizes ) {
	SMF_ASSERT( newnum >= 0 );
	if ( resizes || newnum > m_iSize ) {
		resize( newnum );
	}
	m_iNum = newnum;
}

/*
CList<ItemType>::setGranularity
Sets the base size of the array and resizes the array to match.
*/
/** \brief Define o tamanho da base da matriz e a redimensiona para corresponder.
*  \param granularity: granularidade da lista
*/
template< class ItemType >
SMF_INLINE_FORCED void CList<ItemType>::setGranularity( int newgranularity ) {
	int newsize;

	SMF_ASSERT( newgranularity > 0 );
	m_iGranularity = newgranularity;

	if ( m_pData ) {
		// resize it to the closest level of granularity
		newsize = m_iNum + m_iGranularity - 1;
		newsize -= newsize % m_iGranularity;
		if ( newsize != m_iSize ) {
			resize( newsize );
		}
	}
}

/*
CList<ItemType>::getGranularity
Get the current granularity.
*/
/** \brief retorna a granularidade da lista
* \return A granularidade da lista
*/
template< class ItemType >
SMF_INLINE_FORCED int CList<ItemType>::getGranularity() const {
	return m_iGranularity;
}

/*
CList<ItemType>::condense
Resizes the array to exactly the number of elements it contains or frees up memory if empty.
*/
/** \brief Redimensiona a matriz da lista para o tamanho exato igual ao número de elementos que ela contém
* \note libera memória que estiver livre
*
*/
template< class ItemType >
SMF_INLINE_FORCED void CList<ItemType>::condense() {
	if ( m_pData ) {
		if ( m_iNum ) {
			resize( m_iNum );
		} else {
			clear();
		}
	}
}

/*
CList<ItemType>::resize
Allocates memory for the amount of elements requested while keeping the contents intact.
Contents are copied using their = operator so that data is correnctly instantiated.
*/
/** \brief Aloca memória para a quantidade de elementos requisitada enquanto mantém o conteúdo intacto.
* \note O conteúdo é copiado usando o perador de atribuição "=" do próprio objeto
* \param int newsize: novo tyamanho da lista
*/
template< class ItemType >
SMF_INLINE_FORCED void CList<ItemType>::resize( int newsize ) {
	ItemType	*temp;
	int		i;

	SMF_ASSERT( newsize >= 0 );

	// free up the list if no data is being reserved
	if ( newsize <= 0 ) {
		clear();
		return;
	}

	if ( newsize == m_iSize ) {
		// not changing the size, so just exit
		return;
	}

	temp	= m_pData;
	m_iSize	= newsize;
	if ( m_iSize < m_iNum ) {
		m_iNum = m_iSize;
	}

	// copy the old list into our new one
	m_pData = new ItemType[ m_iSize ];
	for( i = 0; i < m_iNum; i++ ) {
		m_pData[ i ] = temp[ i ];
	}

	// delete the old list if it exists
	if ( temp ) {
		delete[] temp;
	}
}

/*
CList<ItemType>::resize
Allocates memory for the amount of elements requested while keeping the contents intact.
Contents are copied using their = operator so that data is correctly instantiated.
*/
/** \brief Aloca memória para a quantidade de elementos requisitada enquanto mantém o conteúdo intacto.
* \note O conteúdo é copiado usando o perador de atribuição "=" do próprio objeto
* \param int newsize: novo tyamanho da lista
* \param int newgranularity: nova granularidade da lista
*/
template< class ItemType >
SMF_INLINE_FORCED void CList<ItemType>::resize( int newsize, int newgranularity ) {
	ItemType	*temp;
	int		i;

	SMF_ASSERT( newsize >= 0 );

	SMF_ASSERT( newgranularity > 0 );
	m_iGranularity = newgranularity;

	// free up the list if no data is being reserved
	if ( newsize <= 0 ) {
		clear();
		return;
	}

	temp	= m_pData;
	m_iSize	= newsize;
	if ( m_iSize < m_iNum ) {
		m_iNum = m_iSize;
	}

	// copy the old list into our new one
	m_pData = new ItemType[ m_iSize ];
	for( i = 0; i < m_iNum; i++ ) {
		m_pData[ i ] = temp[ i ];
	}

	// delete the old list if it exists
	if ( temp ) {
		delete[] temp;
	}
}

/*
CList<ItemType>::assureSize
Makes sure the list has at least the given number of elements.
*/
/** \brief assegura que a liste tem pelo menos o dado número de elementos
* \param  int newSize: Novo tamanho da lista
*/
template< class ItemType >
SMF_INLINE_FORCED void CList<ItemType>::assureSize( int newSize ) {
	int newNum = newSize;

	if ( newSize > m_iSize ) {

		if ( m_iGranularity == 0 ) {	// this is a hack to fix our memset classes
			m_iGranularity = 16;
		}

		newSize += m_iGranularity - 1;
		newSize -= newSize % m_iGranularity;
		resize( newSize );
	}

	m_iNum = newNum;
}

/*
CList<ItemType>::assureSize
Makes sure the list has at least the given number of elements and initialize any elements not yet initialized.
*/
/** \brief assegura que a lista possui o númerodado e elementos e inicializa-os
* \param  int newSize: novo tamanho da lista
* \param const ItemType &initValue: referência para um objeto que será utilizado como valor inicial
*/
template< class ItemType >
SMF_INLINE_FORCED void CList<ItemType>::assureSize( int newSize, const ItemType &initValue ) {
	int newNum = newSize;

	if ( newSize > m_iSize ) {

		if ( m_iGranularity == 0 ) {	// this is a hack to fix our memset classes
			m_iGranularity = 16;
		}

		newSize += m_iGranularity - 1;
		newSize -= newSize % m_iGranularity;
		m_iNum = m_iSize;
		resize( newSize );

		for ( int i = m_iNum; i < newSize; i++ ) {
			m_pData[i] = initValue;
		}
	}

	m_iNum = newNum;
}

/*
CList<ItemType>::assureSizeAlloc
Makes sure the list has at least the given number of elements and allocates any elements using the allocator.
NOTE: This function can only be called on lists containing pointers. Calling it
on non-pointer lists will cause a compiler error.
*/
/** \brief assegura que a lista tem no mínimo a quantidade de elementos passados
* aloca elementos caso necessário
* \note este método só funciona em listas contendo ponteiros. 
* \note Em listas de objetos que não sejam ponteiros causará erro de compilação
*
*/
template< class ItemType >
SMF_INLINE_FORCED void CList<ItemType>::assureSizeAlloc( int newSize, new_t *allocator ) {
	int newNum = newSize;

	if ( newSize > m_iSize ) {

		if ( m_iGranularity == 0 ) {	// this is a hack to fix our memset classes
			m_iGranularity = 16;
		}

		newSize += m_iGranularity - 1;
		newSize -= newSize % m_iGranularity;
		m_iNum = m_iSize;
		resize( newSize );

		for ( int i = m_iNum; i < newSize; i++ ) {
			m_pData[i] = (*allocator)();
		}
	}

	m_iNum = newNum;
}

/*
CList<ItemType>::operator=
Copies the contents and size attributes of another list.
*/
/** \brief operador de atrinuição "=" ou cópia
* \param  const CList<ItemType> &other: objeto o qual se deseja obter a cópia
*/
template< class ItemType >
SMF_INLINE_FORCED CList<ItemType> &CList<ItemType>::operator=( const CList<ItemType> &other ) {
	int	i;

	clear();

	m_iNum			= other.m_iNum;
	m_iSize		= other.m_iSize;
	m_iGranularity	= other.m_iGranularity;

	if ( m_iSize ) {
		m_pData = new ItemType[ m_iSize ];
		for( i = 0; i < m_iNum; i++ ) {
			m_pData[ i ] = other.m_pData[ i ];
		}
	}

	return *this;
}

/*
CList<ItemType>::operator[] const
Access operator.  Index must be within range or an SMF_ASSERT will be issued in debug builds.
Release builds do no range checking.
*/
/** \brief operador de subscrição "[]"
* \param int index: O índice do objeto na lista
*/
template< class ItemType >
SMF_INLINE_FORCED const ItemType &CList<ItemType>::operator[]( int index ) const {
	SMF_ASSERT( index >= 0 );
	SMF_ASSERT( index < m_iNum );

	return m_pData[ index ];
}

/*
CList<ItemType>::operator[]
Access operator.  Index must be within range or an SMF_ASSERT will be issued in debug builds.
Release builds do no range checking.
*/
/** \brief operador de subscrição "[]"
* \param int index: O índice do objeto na lista
*/
template< class ItemType >
SMF_INLINE_FORCED ItemType &CList<ItemType>::operator[]( int index ) {
	SMF_ASSERT( index >= 0 );
	SMF_ASSERT( index < m_iNum );

	return m_pData[ index ];
}

/*
CList<ItemType>::getPtr
Returns a pointer to the begining of the array.  Useful for iterating through the list in loops.
Note: may return NULL if the list is empty.
FIXME: Create an iterator template for this kind of thing.
*/
/** \brief Retorna um ponteiro para o início da matriz da Lista
* \return Um ponteiro para o início da Lista
* \return Um ponteiro NULL se a lista estiver vazia
* \note Útil para iteragir com a lista através de loops
*/
template< class ItemType >
SMF_INLINE_FORCED ItemType *CList<ItemType>::getPtr() {
	return m_pData;
}

/*
CList<ItemType>::getPtr
Returns a pointer to the begining of the array.  Useful for iterating through the list in loops.
Note: may return NULL if the list is empty.
FIXME: Create an iterator template for this kind of thing.
*/
/** \brief Retorna um ponteiro para o início da matriz da Lista
* \return Um ponteiro para o início da Lista
* \return Um ponteiro NULL se a lista estiver vazia
* \note Útil para iteragir com a lista através de loops
*/
template< class ItemType >
const SMF_INLINE_FORCED ItemType *CList<ItemType>::getPtr() const {
	return m_pData;
}

/*
CList<ItemType>::alloc
Returns a reference to a new data element at the end of the list.
*/
/** \brief Aloca um novo container para objetos no final da lista
* \return O endereço do novo container adicionado
*/
template< class ItemType >
SMF_INLINE_FORCED ItemType &CList<ItemType>::alloc() {
	if ( !m_pData ) {
		resize( m_iGranularity );
	}

	if ( m_iNum == m_iSize ) {
		resize( m_iSize + m_iGranularity );
	}

	return m_pData[ m_iNum++ ];
}

/*
CList<ItemType>::append
Increases the size of the list by one element and copies the supplied data into it.
Returns the index of the new element.
*/
/** \brief Anexa um objeto na lista
* Copia o objeto passado no parâmetro para a lista e aumenta o tamanho da lista em uma unidade
* \param ItemType const & obj: Uma referência para o objeto ao qual se deseja anexar
* \retorna o índice do novo elemento
*/
template< class ItemType >
SMF_INLINE_FORCED int CList<ItemType>::append( ItemType const & obj ) {
	if ( !m_pData ) {
		resize( m_iGranularity );
	}

	if ( m_iNum == m_iSize ) {
		int newsize;

		if ( m_iGranularity == 0 ) {	// this is a hack to fix our memset classes
			m_iGranularity = 16;
		}
		newsize = m_iSize + m_iGranularity;
		resize( newsize - newsize % m_iGranularity );
	}

	m_pData[ m_iNum ] = obj;
	m_iNum++;

	return m_iNum - 1;
}


/*
CList<ItemType>::insert
Increases the size of the list by at leat one element if necessary 
and inserts the supplied data into it.
Returns the index of the new element.
*/
/** \brief insere o objeto passado na lista
* corrige o tamanho da lista
* \param ItemType const & obj: Uma referência ao objeto a ser adicionado na lista
* \param int index: O indice no qual se deseja inserir o objeto 
* \return o Índice ao qual o objeto foi adicionado
*/
template< class ItemType >
SMF_INLINE_FORCED int CList<ItemType>::insert( ItemType const & obj, int index ) {
	if ( !m_pData ) {
		resize( m_iGranularity );
	}

	if ( m_iNum == m_iSize ) {
		int newsize;

		if ( m_iGranularity == 0 ) {	// this is a hack to fix our memset classes
			m_iGranularity = 16;
		}
		newsize = m_iSize + m_iGranularity;
		resize( newsize - newsize % m_iGranularity );
	}

	if ( index < 0 ) {
		index = 0;
	}
	else if ( index > m_iNum ) {
		index = m_iNum;
	}
	for ( int i = m_iNum; i > index; --i ) {
		m_pData[i] = m_pData[i-1];
	}
	m_iNum++;
	m_pData[index] = obj;
	return index;
}

/*
CList<ItemType>::append
adds the other list to this one
Returns the size of the new combined list
*/
/** \brief junta a lista passada no parâmetro à esta lista
* \param const CList<ItemType> &other: referência para a lista a qual se deseja juntar
* \return Retorna o tamanho da nova lista combinada
*/
template< class ItemType >
SMF_INLINE_FORCED int CList<ItemType>::append( const CList<ItemType> &other ) {
	if ( !m_pData ) {
		if ( m_iGranularity == 0 ) {	// this is a hack to fix our memset classes
			m_iGranularity = 16;
		}
		resize( m_iGranularity );
	}

	int n = other.m_iNum();
	for (int i = 0; i < n; i++) {
		append(other[i]);
	}

	return m_iNum();
}

/*
CList<ItemType>::addUnique
Adds the data to the list if it doesn't already exist.  Returns the index of the data in the list.
*/
/** \brief adiciona o objeto à lista se ele íainda não estiver lá.
* \return Retorna o índice do objeto adicionado
* \param ItemType const & obj: Uma referência ao objeto a ser adicionado
*/
template< class ItemType >
SMF_INLINE_FORCED int CList<ItemType>::addUnique( ItemType const & obj ) {
	int index;

	index = findIndex( obj );
	if ( index < 0 ) {
		index = append( obj );
	}

	return index;
}

/*
CList<ItemType>::findIndex
Searches for the specified data in the m_pData and returns it's index.  Returns -1 if the data is not found.
*/
/** \brief Procura por um objeto específico e retorna seu índice
* \param ItemType const & obj: Uma referência para o objeto ao qual se deseja obter o índice
* \return -1 caso o objeto não seja encontrado
*/
template< class ItemType >
SMF_INLINE_FORCED int CList<ItemType>::findIndex( ItemType const & obj ) const {
	int i;

	for( i = 0; i < m_iNum; i++ ) {
		if ( m_pData[ i ] == obj ) {
			return i;
		}
	}

	// Not found
	return -1;
}

/*
CList<ItemType>::find
Searches for the specified data in the m_pData and returns it's address. Returns NULL if the data is not found.
*/
/** \brief procura um objeto na lista
* \param  ItemType const & obj: Uma referência pdo objeto que se deseja encontrar na lista
* \return Um ponteiro para o elemento da lista que contém o objeto ou NULL caso não encontre o objeto
*/
template< class ItemType >
SMF_INLINE_FORCED ItemType *CList<ItemType>::find( ItemType const & obj ) const {
	int i;

	i = findIndex( obj );
	if ( i >= 0 ) {
		return &m_pData[ i ];
	}

	return NULL;
}

/*
CList<ItemType>::findNull
Searches for a NULL pointer in the list.  Returns -1 if NULL is not found.
NOTE: This function can only be called on lists containing pointers. Calling it
on non-pointer lists will cause a compiler error.
*/
/** \brief Procura um ponteiro nulo na lista
* \return -1 caso o ponteiro nulo não seja encontrado
* \note Este método só pode ser invocado em lista que contém ponteiros.
* se o método for invocado numa lista que não contém ponteiros  ocorrerá erro de compilação
*/
template< class ItemType >
SMF_INLINE_FORCED int CList<ItemType>::findNull() const {
	int i;

	for( i = 0; i < m_iNum; i++ ) {
		if ( m_pData[ i ] == NULL ) {
			return i;
		}
	}

	// Not found
	return -1;
}

/*
CList<ItemType>::getIndexOf
Takes a pointer to an element in the list and returns the index of the element.
This is NOT a guarantee that the object is really in the list. 
Function will SMF_ASSERT in debug builds if pointer is outside the bounds of the list,
but remains silent in release builds.
*/
/** \brief retorna o índice do elemento 
* Recebe como parâmetro um ponteiro para um objeto da lista e retorna seu índice.
* \param ItemType const *objptr: Um ponteiro para o objeto que se deseja saber o índice
* \return O valor do índice, ou -1 caso o objeto não esteja na lista
* \note Não há garantia de que o objeto esteja realmente na lista.
*/
template< class ItemType >
SMF_INLINE_FORCED int CList<ItemType>::getIndexOf( ItemType const *objptr ) const {
	int index;

	index = objptr - m_pData;
	if ( index <= 0 || index > m_iNum) return -1;

	return index;
}

/*
CList<ItemType>::removeIndex
Removes the element at the specified index and moves all data following the element down to fill in the gap.
The number of elements in the list is reduced by one.  Returns false if the index is outside the bounds of the list.
Note that the element is not destroyed, so any memory used by it may not be freed until the destruction of the list.
*/
/** \brief remove um objeto do ´pindice especificado no parâmetro passado
* remove o objeto se ele estiver na lista e move todos os objetos que o seguem para não deixar buraco.
* O número de elementos da lista é reduzido em uma unidade.
* \return falso se o índice estiver fora dos limites da lista
* \note O objeto não é destruído, portanto a memória utilizada por ele não é liberada até a destruição da lista.
*/
template< class ItemType >
SMF_INLINE_FORCED bool CList<ItemType>::removeIndex( int index ) {
	int i;

	SMF_ASSERT( m_pData != NULL );
	SMF_ASSERT( index >= 0 );
	SMF_ASSERT( index < m_iNum );

	if ( ( index < 0 ) || ( index >= m_iNum ) ) {
		return false;
	}

	m_iNum--;
	for( i = index; i < m_iNum; i++ ) {
		m_pData[ i ] = m_pData[ i + 1 ];
	}

	return true;
}

/*

CList<ItemType>::remove
Removes the element if it is found within the list and moves all data following the element down to fill in the gap.
The number of elements in the list is reduced by one.  Returns false if the data is not found in the list.  Note that
the element is not destroyed, so any memory used by it may not be freed until the destruction of the list.
*/
/** \brief remove um objeto da lista
* \param  ItemType const & obj: uma referência para do objeto a ser removido
* remove o objeto se ele estiver na lista e move todos os objetos que o seguem para não deixar buraco.
* O número de elementos da lista é reduzido em uma unidade.
* \return falso se o objeto não foi encontrado na lista
* \note O objeto não é destruído, portanto a memória utilizada por ele não é liberada até a destruição da lista.
*/
template< class ItemType >
SMF_INLINE_FORCED bool CList<ItemType>::remove( ItemType const & obj ) {
	int index;

	index = findIndex( obj );
	if ( index >= 0 ) {
		return removeIndex( index );
	}
	
	return false;
}

/*
CList<ItemType>::sort
Performs a qsort on the list using the supplied comparison function.  Note that the data is merely moved around the
list, so any pointers to data within the list may no longer be valid.
*/
/** \brief ordena a lista utilizando a função de comparação passada no parâmetro
* \param cmp_t *compare : ponteiro para a função de comparação
* \note Os dados são meramente movidos dentro da lista, portanto, vários ponteiros para dado podem não ser válidos
*/
template< class ItemType >
SMF_INLINE_FORCED void CList<ItemType>::sort( cmp_t *compare ) {
	if ( !m_pData ) {
		return;
	}
	typedef int cmp_c(const void *, const void *);

	cmp_c *vCompare = (cmp_c *)compare;
	qsort( ( void * )m_pData, ( size_t )m_iNum, sizeof( ItemType ), vCompare );
}

/*
CList<ItemType>::sortSubSection
Sorts a subsection of the list.
*/
/** \brief Ordena uma parte da lista
* \param  int startIndex: a posicão inicial da lista a qual se deseja iniciar a ordenação
* \param  int endIndex: a posição final (índice) da lista , na qual a oredenação deve acabar
* \param  cmp_t *compare: a função de comparação dos dois objetos, a qual decidirá se a ordem de um objeto é maior ou menor que do outro
*/
template< class ItemType >
SMF_INLINE_FORCED void CList<ItemType>::sortSubSection( int startIndex, int endIndex, cmp_t *compare ) {
	if ( !m_pData ) {
		return;
	}
	if ( startIndex < 0 ) {
		startIndex = 0;
	}
	if ( endIndex >= m_iNum ) {
		endIndex = m_iNum - 1;
	}
	if ( startIndex >= endIndex ) {
		return;
	}
	typedef int cmp_c(const void *, const void *);

	cmp_c *vCompare = (cmp_c *)compare;
	qsort( ( void * )( &m_pData[startIndex] ), ( size_t )( endIndex - startIndex + 1 ), sizeof( ItemType ), vCompare );
}

/*
CList<ItemType>::swapList
Swaps the contents of two lists
*/
/** \brief Troca o conteúde de duas listas entre si
* \param other: A lista a qual se deseja realizar a troca
*/
template< class ItemType >
SMF_INLINE_FORCED void CList<ItemType>::swapList( CList<ItemType> &other ) {
	swapElements( m_iNum, other.m_iNum );
	swapElements( m_iSize, other.m_iSize );
	swapElements( m_iGranularity, other.m_iGranularity );
	swapElements( m_pData, other.m_pData );
}

} //end SMF
#endif /* !__SMF_LIST_H__*/
