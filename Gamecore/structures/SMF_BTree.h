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


#ifndef _SMF__BTREE_H__
#define _SMF__BTREE_H__
#include "../SMF_Config.h"
#include "../util/SMF_HeapBlockAloc.h"
#include "../util/SMF_UtilStructs.h"
#include "../exceptions/all.h"
namespace SMF{
using namespace Util;


//#define BTREE_CHECK
/**
 * \class CBtreeBalCBTreeNode
 *
 * \ingroup SMF_Data_Structures
 *
 * \if pt_br
 * \brief Nó da Classe \ref CBtree \ref
 * \elseif us_en
 * \brief Node for the class \ref CBtree \ref
 * \endif
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
template< class objType, class keyType >
class SMF_API CBTreeNode {
public:
	keyType							key;			// key used for sorting
	objType *						object;			// if != NULL pointer to object stored in leaf node
	CBTreeNode *					parent;			// parent node
	CBTreeNode *					next;			// next sibling
	CBTreeNode *					prev;			// prev sibling
	int								numChildren;	// number of children
	CBTreeNode *					firstChild;		// first child
	CBTreeNode *					lastChild;		// last child
};

/**
 * \class CBtreeBal
 *
 * \ingroup SMF_Data_Structures
 *
 * \if pt_br
 * \brief Implementa uma Árvore Binária de Busca Balanceada
 * \elseif us_en
 * \brief Balanced Search Tree
 * \endif
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \see http://en.wikipedia.org/wiki/Self-balancing_binary_search_tree
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
template< class objType, class keyType, int maxChildrenPerNode >
class SMF_API CBtreeBal {
public:
									CBtreeBal( void );
									~CBtreeBal( void );

	void							init( void );
	void							shutdown( void );
	/// add an object to the tree
	CBTreeNode<objType,keyType> *	add( objType *object, keyType key );						
	/// remove an object node from the tree
	void							remove( CBTreeNode<objType,keyType> *node );				
	/// find an object using the given key
	objType *						find( keyType key ) const;									
	/// find an object with the smallest key larger equal the given key
	objType *						findSmallestLargerEqual( keyType key ) const;	
	/// find an object with the largest key smaller equal the given key
	objType *						findLargestSmallerEqual( keyType key ) const;				
	/// returns the root node of the tree
	CBTreeNode<objType,keyType> *	getRoot( void ) const;
	/// returns the total number of nodes in the tree
	int								getNodeCount( void ) const;									
	/// goes through all nodes of the tree
	CBTreeNode<objType,keyType> *	getNext( CBTreeNode<objType,keyType> *node ) const;		
	/// goes through all leaf nodes of the tree
	CBTreeNode<objType,keyType> *	getNextLeaf( CBTreeNode<objType,keyType> *node ) const;	

private:
	CBTreeNode<objType,keyType> *	root;
	CBlockAlloc<CBTreeNode<objType,keyType>,128>	nodeAllocator;

	CBTreeNode<objType,keyType> *	AllocNode( void );
	void							FreeNode( CBTreeNode<objType,keyType> *node );
	void							SplitNode( CBTreeNode<objType,keyType> *node );
	CBTreeNode<objType,keyType> *	MergeNodes( CBTreeNode<objType,keyType> *node1, CBTreeNode<objType,keyType> *node2 );

	void							CheckTree_r( CBTreeNode<objType,keyType> *node, int &numNodes ) const;
	void							CheckTree( void ) const;
};

#if !defined(_MSC_VER) && !defined(SMF_INLINE_FORCED)
 #define SMF_INLINE_FORCED inline
#endif

template< class objType, class keyType, int maxChildrenPerNode >
SMF_INLINE_FORCED CBtreeBal<objType,keyType,maxChildrenPerNode>::CBtreeBal( void ) {
	SMF_ASSERT( maxChildrenPerNode >= 4 );
	root = NULL;
}

template< class objType, class keyType, int maxChildrenPerNode >
SMF_INLINE_FORCED CBtreeBal<objType,keyType,maxChildrenPerNode>::~CBtreeBal( void ) {
	shutdown();
}

template< class objType, class keyType, int maxChildrenPerNode >
SMF_INLINE_FORCED void CBtreeBal<objType,keyType,maxChildrenPerNode>::init( void ) {
	root = AllocNode();
}

template< class objType, class keyType, int maxChildrenPerNode >
SMF_INLINE_FORCED void CBtreeBal<objType,keyType,maxChildrenPerNode>::shutdown( void ) {
	nodeAllocator.shutdown();
	root = NULL;
}

template< class objType, class keyType, int maxChildrenPerNode >
SMF_INLINE_FORCED CBTreeNode<objType,keyType> *CBtreeBal<objType,keyType,maxChildrenPerNode>::add( objType *object, keyType key ) {
	CBTreeNode<objType,keyType> *node, *child, *newNode;

	if ( root->numChildren >= maxChildrenPerNode ) {
		newNode = AllocNode();
		newNode->key = root->key;
		newNode->firstChild = root;
		newNode->lastChild = root;
		newNode->numChildren = 1;
		root->parent = newNode;
		SplitNode( root );
		root = newNode;
	}

	newNode = AllocNode();
	newNode->key = key;
	newNode->object = object;

	for ( node = root; node->firstChild != NULL; node = child ) {

		if ( key > node->key ) {
			node->key = key;
		}

		// find the first child with a key larger equal to the key of the new node
		for( child = node->firstChild; child->next; child = child->next ) {
			if ( key <= child->key ) {
				break;
			}
		}

		if ( child->object ) {

			if ( key <= child->key ) {
				// insert new node before child
				if ( child->prev ) {
					child->prev->next = newNode;
				} else {
					node->firstChild = newNode;
				}
				newNode->prev = child->prev;
				newNode->next = child;
				child->prev = newNode;
			} else {
				// insert new node after child
				if ( child->next ) {
					child->next->prev = newNode;
				} else {
					node->lastChild = newNode;
				}
				newNode->prev = child;
				newNode->next = child->next;
				child->next = newNode;
			}

			newNode->parent = node;
			node->numChildren++;

#ifdef BTREE_CHECK
			CheckTree();
#endif

			return newNode;
		}

		// make sure the child has room to store another node
		if ( child->numChildren >= maxChildrenPerNode ) {
			SplitNode( child );
			if ( key <= child->prev->key ) {
				child = child->prev;
			}
		}
	}

	// we only end up here if the root node is empty
	newNode->parent = root;
	root->key = key;
	root->firstChild = newNode;
	root->lastChild = newNode;
	root->numChildren++;

#ifdef BTREE_CHECK
	CheckTree();
#endif

	return newNode;
}

template< class objType, class keyType, int maxChildrenPerNode >
SMF_INLINE_FORCED void CBtreeBal<objType,keyType,maxChildrenPerNode>::remove( CBTreeNode<objType,keyType> *node ) {
	CBTreeNode<objType,keyType> *parent;

	SMF_ASSERT( node->object != NULL );

	// unlink the node from it's parent
	if ( node->prev ) {
		node->prev->next = node->next;
	} else {
		node->parent->firstChild = node->next;
	}
	if ( node->next ) {
		node->next->prev = node->prev;
	} else {
		node->parent->lastChild = node->prev;
	}
	node->parent->numChildren--;

	// make sure there are no parent nodes with a single child
	for ( parent = node->parent; parent != root && parent->numChildren <= 1; parent = parent->parent ) {

		if ( parent->next ) {
			parent = MergeNodes( parent, parent->next );
		} else if ( parent->prev ) {
			parent = MergeNodes( parent->prev, parent );
		}

		// a parent may not use a key higher than the key of it's last child
		if ( parent->key > parent->lastChild->key ) {
			parent->key = parent->lastChild->key;
		}

		if ( parent->numChildren > maxChildrenPerNode ) {
			SplitNode( parent );
			break;
		}
	}
	for ( ; parent != NULL && parent->lastChild != NULL; parent = parent->parent ) {
		// a parent may not use a key higher than the key of it's last child
		if ( parent->key > parent->lastChild->key ) {
			parent->key = parent->lastChild->key;
		}
	}

	// free the node
	FreeNode( node );

	// remove the root node if it has a single internal node as child
	if ( root->numChildren == 1 && root->firstChild->object == NULL ) {
		CBTreeNode<objType,keyType> *oldRoot = root;
		root->firstChild->parent = NULL;
		root = root->firstChild;
		FreeNode( oldRoot );
	}

#ifdef BTREE_CHECK
	CheckTree();
#endif
}

template< class objType, class keyType, int maxChildrenPerNode >
SMF_INLINE_FORCED objType *CBtreeBal<objType,keyType,maxChildrenPerNode>::find( keyType key ) const {
	CBTreeNode<objType,keyType> *node;

	for ( node = root->firstChild; node != NULL; node = node->firstChild ) {
		while( node->next ) {
			if ( node->key >= key ) {
				break;
			}
			node = node->next;
		}
		if ( node->object ) {
			if ( node->key == key ) {
				return node->object;
			} else {
				return NULL;
			}
		}
	}
	return NULL;
}

template< class objType, class keyType, int maxChildrenPerNode >
SMF_INLINE_FORCED objType *CBtreeBal<objType,keyType,maxChildrenPerNode>::findSmallestLargerEqual( keyType key ) const {
	CBTreeNode<objType,keyType> *node;

	for ( node = root->firstChild; node != NULL; node = node->firstChild ) {
		while( node->next ) {
			if ( node->key >= key ) {
				break;
			}
			node = node->next;
		}
		if ( node->object ) {
			if ( node->key >= key ) {
				return node->object;
			} else {
				return NULL;
			}
		}
	}
	return NULL;
}

template< class objType, class keyType, int maxChildrenPerNode >
SMF_INLINE_FORCED objType *CBtreeBal<objType,keyType,maxChildrenPerNode>::findLargestSmallerEqual( keyType key ) const {
	CBTreeNode<objType,keyType> *node;

	for ( node = root->lastChild; node != NULL; node = node->lastChild ) {
		while( node->prev ) {
			if ( node->key <= key ) {
				break;
			}
			node = node->prev;
		}
		if ( node->object ) {
			if ( node->key <= key ) {
				return node->object;
			} else {
				return NULL;
			}
		}
	}
	return NULL;
}

template< class objType, class keyType, int maxChildrenPerNode >
SMF_INLINE_FORCED CBTreeNode<objType,keyType> *CBtreeBal<objType,keyType,maxChildrenPerNode>::getRoot( void ) const {
	return root;
}

template< class objType, class keyType, int maxChildrenPerNode >
SMF_INLINE_FORCED int CBtreeBal<objType,keyType,maxChildrenPerNode>::getNodeCount( void ) const {
	return nodeAllocator.getAllocCount();
}

template< class objType, class keyType, int maxChildrenPerNode >
SMF_INLINE_FORCED CBTreeNode<objType,keyType> *CBtreeBal<objType,keyType,maxChildrenPerNode>::getNext( CBTreeNode<objType,keyType> *node ) const {
	if ( node->firstChild ) {
		return node->firstChild;
	} else {
		while( node && node->next == NULL ) {
			node = node->parent;
		}
		return node;
	}
}

template< class objType, class keyType, int maxChildrenPerNode >
SMF_INLINE_FORCED CBTreeNode<objType,keyType> *CBtreeBal<objType,keyType,maxChildrenPerNode>::getNextLeaf( CBTreeNode<objType,keyType> *node ) const {
	if ( node->firstChild ) {
		while ( node->firstChild ) {
			node = node->firstChild;
		}
		return node;
	} else {
		while( node && node->next == NULL ) {
			node = node->parent;
		}
		if ( node ) {
			node = node->next;
			while ( node->firstChild ) {
				node = node->firstChild;
			}
			return node;
		} else {
			return NULL;
		}
	}
}

template< class objType, class keyType, int maxChildrenPerNode >
SMF_INLINE_FORCED CBTreeNode<objType,keyType> *CBtreeBal<objType,keyType,maxChildrenPerNode>::AllocNode( void ) {
	CBTreeNode<objType,keyType> *node = nodeAllocator.alloc();
	node->key = 0;
	node->parent = NULL;
	node->next = NULL;
	node->prev = NULL;
	node->numChildren = 0;
	node->firstChild = NULL;
	node->lastChild = NULL;
	node->object = NULL;
	return node;
}

template< class objType, class keyType, int maxChildrenPerNode >
SMF_INLINE_FORCED void CBtreeBal<objType,keyType,maxChildrenPerNode>::FreeNode( CBTreeNode<objType,keyType> *node ) {
	nodeAllocator.free( node );
}

template< class objType, class keyType, int maxChildrenPerNode >
SMF_INLINE_FORCED void CBtreeBal<objType,keyType,maxChildrenPerNode>::SplitNode( CBTreeNode<objType,keyType> *node ) {
	int i;
	CBTreeNode<objType,keyType> *child, *newNode;

	// allocate a new node
	newNode = AllocNode();
	newNode->parent = node->parent;

	// divide the children over the two nodes
	child = node->firstChild;
	child->parent = newNode;
	for ( i = 3; i < node->numChildren; i += 2 ) {
		child = child->next;
		child->parent = newNode;
	}

	newNode->key = child->key;
	newNode->numChildren = node->numChildren / 2;
	newNode->firstChild = node->firstChild;
	newNode->lastChild = child;

	node->numChildren -= newNode->numChildren;
	node->firstChild = child->next;

	child->next->prev = NULL;
	child->next = NULL;

	// add the new child to the parent before the split node
	SMF_ASSERT( node->parent->numChildren < maxChildrenPerNode );

	if ( node->prev ) {
		node->prev->next = newNode;
	} else {
		node->parent->firstChild = newNode;
	}
	newNode->prev = node->prev;
	newNode->next = node;
	node->prev = newNode;

	node->parent->numChildren++;
}

template< class objType, class keyType, int maxChildrenPerNode >
SMF_INLINE_FORCED CBTreeNode<objType,keyType> *CBtreeBal<objType,keyType,maxChildrenPerNode>::MergeNodes( CBTreeNode<objType,keyType> *node1, CBTreeNode<objType,keyType> *node2 ) {
	CBTreeNode<objType,keyType> *child;

	SMF_ASSERT( node1->parent == node2->parent );
	SMF_ASSERT( node1->next == node2 && node2->prev == node1 );
	SMF_ASSERT( node1->object == NULL && node2->object == NULL );
	SMF_ASSERT( node1->numChildren >= 1 && node2->numChildren >= 1 );

	for ( child = node1->firstChild; child->next; child = child->next ) {
		child->parent = node2;
	}
	child->parent = node2;
	child->next = node2->firstChild;
	node2->firstChild->prev = child;
	node2->firstChild = node1->firstChild;
	node2->numChildren += node1->numChildren;

	// unlink the first node from the parent
	if ( node1->prev ) {
		node1->prev->next = node2;
	} else {
		node1->parent->firstChild = node2;
	}
	node2->prev = node1->prev;
	node2->parent->numChildren--;

	FreeNode( node1 );

	return node2;
}

template< class objType, class keyType, int maxChildrenPerNode >
SMF_INLINE_FORCED void CBtreeBal<objType,keyType,maxChildrenPerNode>::CheckTree_r( CBTreeNode<objType,keyType> *node, int &numNodes ) const {
	int numChildren;
	CBTreeNode<objType,keyType> *child;

	numNodes++;

	// the root node may have zero children and leaf nodes always have zero children, all other nodes should have at least 2 and at most maxChildrenPerNode children
	SMF_ASSERT( ( node == root ) || ( node->object != NULL && node->numChildren == 0 ) || ( node->numChildren >= 2 && node->numChildren <= maxChildrenPerNode ) );
	// the key of a node may never be larger than the key of it's last child
	SMF_ASSERT( ( node->lastChild == NULL ) || ( node->key <= node->lastChild->key ) );

	numChildren = 0;
	for ( child = node->firstChild; child; child = child->next ) {
		numChildren++;
		// make sure the children are properly linked
		if ( child->prev == NULL ) {
			SMF_ASSERT( node->firstChild == child );
		} else {
			SMF_ASSERT( child->prev->next == child );
		}
		if ( child->next == NULL ) {
			SMF_ASSERT( node->lastChild == child );
		} else {
			SMF_ASSERT( child->next->prev == child );
		}
		// recurse down the tree
		CheckTree_r( child, numNodes );
	}
	// the number of children should equal the number of linked children
	SMF_ASSERT( numChildren == node->numChildren );
}

template< class objType, class keyType, int maxChildrenPerNode >
SMF_INLINE_FORCED void CBtreeBal<objType,keyType,maxChildrenPerNode>::CheckTree( void ) const {
	int numNodes = 0;
	CBTreeNode<objType,keyType> *node, *lastNode;

	CheckTree_r( root, numNodes );

	// the number of nodes in the tree should equal the number of allocated nodes
	SMF_ASSERT( numNodes == nodeAllocator.getAllocCount() );

	// all the leaf nodes should be ordered
	lastNode = getNextLeaf( getRoot() );
	if ( lastNode ) {
		for ( node = getNextLeaf( lastNode ); node; lastNode = node, node = getNextLeaf( node ) ) {
			SMF_ASSERT( lastNode->key <= node->key );
		}
	}
}

} //end SMF
#endif /* !__BTREE_H__ */
