//STARTHEADER
// $Id: SearchTree.hh 3107 2013-05-03 15:47:47Z salam $
//
// Copyright (c) 2005-2011, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of FastJet.
//
//  FastJet is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  The algorithms that underlie FastJet have required considerable
//  development and are described in hep-ph/0512210. If you use
//  FastJet as part of work towards a scientific publication, please
//  include a citation to the FastJet paper.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//ENDHEADER


#ifndef __FASTJET_SEARCHTREE_HH__
#define __FASTJET_SEARCHTREE_HH__

#include<vector>
#include<cassert>
#include<cstddef>
#include "fastjet/internal/base.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


//======================================================================
/// \if internal_doc
/// @ingroup internal
/// \class SearchTree
/// Efficient class for a search tree
///
/// This is the class for a search tree designed to be especially efficient
/// when looking for successors and predecessors (to be used in Chan's
/// CP algorithm). It has the requirement that the maximum size of the
/// search tree must be known in advance.
/// \endif
template<class T> class SearchTree {
public:

  class Node;
  class circulator;
  class const_circulator;

  /// constructor for a search tree from an ordered vector
  SearchTree(const std::vector<T> & init);

  /// constructor for a search tree from an ordered vector allowing
  /// for future growth beyond the current size, up to max_size
  SearchTree(const std::vector<T> & init, unsigned int max_size);

  /// remove the node corresponding to node_index from the search tree
  void remove(unsigned node_index);
  void remove(typename SearchTree::Node * node);
  void remove(typename SearchTree::circulator & circ);

  /// insert the supplied value into the tree and return a pointer to
  /// the relevant SearchTreeNode.
  //Node * insert(const T & value);
  circulator insert(const T & value);

  const Node & operator[](int i) const {return _nodes[i];};

  /// return the number of elements currently in the search tree
  unsigned int size() const {return _nodes.size() - _available_nodes.size();}

  /// check that the structure we've obtained makes sense...
  void verify_structure();
  void verify_structure_linear() const;
  void verify_structure_recursive(const Node * , const Node * , const Node * ) const;

  /// print out all elements...
  void print_elements();

  // tracking the depth may have some speed overhead -- so leave it 
  // out for the time being...
#ifdef __FASTJET_SEARCHTREE_TRACK_DEPTH
  /// the max depth the tree has ever reached
  inline unsigned int max_depth() const {return _max_depth;};
#else
  inline unsigned int max_depth() const {return 0;};
#endif

  int loc(const Node * node) const ;

  /// return predecessor by walking through the tree
  Node * _find_predecessor(const Node *);
  /// return successor by walking through the tree
  Node * _find_successor(const Node *);

  const Node & operator[](unsigned int i) const {return _nodes[i];};

  /// return a circulator to some place in the tree (with a circulator
  /// you don't care where...)
  const_circulator somewhere() const;
  circulator somewhere();

private:
  
  void _initialize(const std::vector<T> & init);

  std::vector<Node> _nodes;
  std::vector<Node *> _available_nodes;
  Node * _top_node;
  unsigned int _n_removes;

  
  /// recursive routine for doing the initial connections assuming things
  /// are ordered. Assumes this_one's parent is labelled, and was
  /// generated at a scale "scale" -- connections will be carried out
  /// including left edge and excluding right edge
  void _do_initial_connections(unsigned int this_one, unsigned int scale,
			       unsigned int left_edge, unsigned int right_edge,
			       unsigned int depth);

  
#ifdef __FASTJET_SEARCHTREE_TRACK_DEPTH
  unsigned int _max_depth;
#endif

};


//======================================================================
/// \if internal_doc
/// @ingroup internal
/// \class SearchTree::Node
/// A node in the search tree
/// \endif
template<class T> class SearchTree<T>::Node{
public:
  Node() {}; /// default constructor
  
  
  /// returns tree if all the tree-related links are set to null for this node
  bool treelinks_null() const {
    return ((parent==0) && (left==0) && (right==0));};
  
  /// set all the tree-related links are set to null for this node
  inline void nullify_treelinks() {
    parent = NULL; 
    left   = NULL; 
    right  = NULL;
  };
  
  /// if my parent exists, determine whether I am it's left or right
  /// node and set the relevant link equal to XX.
  void reset_parents_link_to_me(Node * XX);
  
  T      value;
  Node * left;
  Node * right;
  Node * parent;
  Node * successor;
  Node * predecessor;
};

//----------------------------------------------------------------------
template<class T> void SearchTree<T>::Node::reset_parents_link_to_me(typename SearchTree<T>::Node * XX) {
  if (parent == NULL) {return;}
  if (parent->right == this) {parent->right = XX;}
  else {parent->left = XX;}
}



//======================================================================
/// \if internal_doc
/// @ingroup internal
/// \class SearchTree::circulator
/// circulator for the search tree
/// \endif
template<class T> class SearchTree<T>::circulator{
public:

  // so that it can access out _node object;
  // note: "class U" needed for clang (v1.1 branches/release_27) compilation
  template<class U> friend class SearchTree<U>::const_circulator;
  friend class SearchTree<T>;

  circulator() : _node(NULL) {}

  circulator(Node * node) : _node(node) {}

  const T * operator->() const {return &(_node->value);}
  T * operator->() {return &(_node->value);}
  const T & operator*() const {return _node->value;}
  T & operator*() {return _node->value;}

  /// prefix increment (structure copied from stl_bvector.h)
  circulator & operator++() {
    _node = _node->successor; 
    return *this;}

  /// postfix increment ["int" argument tells compiler it's postfix]
  /// (structure copied from stl_bvector.h)
  circulator operator++(int) {
    circulator tmp = *this;
    _node = _node->successor; 
    return tmp;}

  /// prefix decrement (structure copied from stl_bvector.h)
  circulator & operator--() {
    _node = _node->predecessor; 
    return *this;}

  /// postfix decrement ["int" argument tells compiler it's postfix]
  /// (structure copied from stl_bvector.h)
  circulator operator--(int) {
    circulator tmp = *this;
    _node = _node->predecessor; 
    return tmp;}

  /// return a circulator referring to the next node
  circulator next() const {
    return circulator(_node->successor);}

  /// return a circulator referring to the previous node
  circulator previous() const {
    return circulator(_node->predecessor);}

  bool operator!=(const circulator & other) const {return other._node != _node;}
  bool operator==(const circulator & other) const {return other._node == _node;}

private:
  Node * _node;
};


//======================================================================
/// \if internal_doc
/// @ingroup internal
/// \class SearchTree::const_circulator
/// A const_circulator for the search tree
/// \endif
template<class T> class SearchTree<T>::const_circulator{
public:

  const_circulator() : _node(NULL) {}

  const_circulator(const Node * node) : _node(node) {}
  const_circulator(const circulator & circ) :_node(circ._node) {}

  const T * operator->() {return &(_node->value);}
  const T & operator*() const {return _node->value;}

  /// prefix increment (structure copied from stl_bvector.h)
  const_circulator & operator++() {
    _node = _node->successor; 
    return *this;}

  /// postfix increment ["int" argument tells compiler it's postfix]
  /// (structure copied from stl_bvector.h)
  const_circulator operator++(int) {
    const_circulator tmp = *this;
    _node = _node->successor; 
    return tmp;}


  /// prefix decrement (structure copied from stl_bvector.h)
  const_circulator & operator--() {
    _node = _node->predecessor; 
    return *this;}

  /// postfix decrement ["int" argument tells compiler it's postfix]
  /// (structure copied from stl_bvector.h)
  const_circulator operator--(int) {
    const_circulator tmp = *this;
    _node = _node->predecessor; 
    return tmp;}

  /// return a circulator referring to the next node
  const_circulator next() const {
    return const_circulator(_node->successor);}

  /// return a circulator referring to the previous node
  const_circulator previous() const {
    return const_circulator(_node->predecessor);}



  bool operator!=(const const_circulator & other) const {return other._node != _node;}
  bool operator==(const const_circulator & other) const {return other._node == _node;}

private:
  const Node * _node;
};




//----------------------------------------------------------------------
/// initialise from a sorted initial array allowing for a larger
/// maximum size of the array...
template<class T> SearchTree<T>::SearchTree(const std::vector<T> & init,
					    unsigned int max_size) :
  _nodes(max_size) {

  _available_nodes.reserve(max_size);
  _available_nodes.resize(max_size - init.size());
  for (unsigned int i = init.size(); i < max_size; i++) {
    _available_nodes[i-init.size()] = &(_nodes[i]);
  }

  _initialize(init);
}

//----------------------------------------------------------------------
/// initialise from a sorted initial array
template<class T> SearchTree<T>::SearchTree(const std::vector<T> & init) :
  _nodes(init.size()), _available_nodes(0) {

  // reserve space for the list of available nodes
  _available_nodes.reserve(init.size());
  _initialize(init);
}

//----------------------------------------------------------------------
/// do the actual hard work of initialization
template<class T> void SearchTree<T>::_initialize(const std::vector<T> & init) {

  _n_removes = 0;
  unsigned n = init.size();
  assert(n>=1);

  // reserve space for the list of available nodes
  //_available_nodes.reserve();

#ifdef __FASTJET_SEARCHTREE_TRACK_DEPTH
  _max_depth     = 0;
#endif


  // validate the input
  for (unsigned int i = 1; i<n; i++) {
    assert(!(init[i] < init[i-1]));
  }
  
  // now initialise the vector; link neighbours in the sequence
  for(unsigned int i = 0; i < n; i++) {
    _nodes[i].value = init[i];
    _nodes[i].predecessor = (& (_nodes[i])) - 1;
    _nodes[i].successor   = (& (_nodes[i])) + 1;
    _nodes[i].nullify_treelinks();
  }
  // make a loop structure so that we can circulate...
  _nodes[0].predecessor = (& (_nodes[n-1]));
  _nodes[n-1].successor = (& (_nodes[0]));

  // now label the rest of the nodes
  unsigned int scale = (n+1)/2;
  unsigned int top   = std::min(n-1,scale);
  _nodes[top].parent = NULL;
  _top_node = &(_nodes[top]);
  _do_initial_connections(top, scale, 0, n, 0);

  // make sure things are sensible...
  //verify_structure();
}



//----------------------------------------------------------------------
template<class T> inline  int SearchTree<T>::loc(const Node * node) const {return node == NULL? 
      -999 : node - &(_nodes[0]);}


//----------------------------------------------------------------------
/// Recursive creation of connections, assuming the _nodes vector is
/// completely filled and ordered
template<class T> void SearchTree<T>::_do_initial_connections(
                                         unsigned int this_one, 
					 unsigned int scale,
					 unsigned int left_edge,
					 unsigned int right_edge,
					 unsigned int depth
					 ) {

#ifdef __FASTJET_SEARCHTREE_TRACK_DEPTH
  // keep track of tree depth for checking things stay reasonable...
  _max_depth = max(depth, _max_depth);
#endif

  //std::cout << this_one << " "<< scale<< std::endl;
  unsigned int ref_new_scale = (scale+1)/2;

  // work through children to our left
  unsigned new_scale = ref_new_scale;
  bool     did_child  = false;
  while(true) {
    int left = this_one - new_scale; // be careful here to use signed int...
    // if there is something unitialised to our left, link to it
    if (left >= static_cast<int>(left_edge) 
	                && _nodes[left].treelinks_null() ) {
      _nodes[left].parent = &(_nodes[this_one]);
      _nodes[this_one].left = &(_nodes[left]);
      // create connections between left_edge and this_one
      _do_initial_connections(left, new_scale, left_edge, this_one, depth+1);
      did_child = true;
      break;
    }
    // reduce the scale so as to try again
    unsigned int old_new_scale = new_scale;
    new_scale = (old_new_scale + 1)/2;
    // unless we've reached end of tree
    if (new_scale == old_new_scale) break;
  }
  if (!did_child) {_nodes[this_one].left = NULL;}


  // work through children to our right
  new_scale = ref_new_scale;
  did_child  = false;
  while(true) {
    unsigned int right = this_one + new_scale;
    if (right < right_edge  && _nodes[right].treelinks_null()) {
      _nodes[right].parent = &(_nodes[this_one]);
      _nodes[this_one].right = &(_nodes[right]);
      // create connections between this_one+1 and right_edge
      _do_initial_connections(right, new_scale, this_one+1,right_edge,depth+1);
      did_child = true;
      break;
    }
    // reduce the scale so as to try again
    unsigned int old_new_scale = new_scale;
    new_scale = (old_new_scale + 1)/2;
    // unless we've reached end of tree
    if (new_scale == old_new_scale) break;
  }
  if (!did_child) {_nodes[this_one].right = NULL;}

}



//----------------------------------------------------------------------
template<class T> void SearchTree<T>::remove(unsigned int node_index) {
  remove(&(_nodes[node_index]));
}

//----------------------------------------------------------------------
template<class T> void SearchTree<T>::remove(circulator & circ) {
  remove(circ._node);
}

//----------------------------------------------------------------------
// Useful reference for this:
//   http://en.wikipedia.org/wiki/Binary_search_tree#Deletion
template<class T> void SearchTree<T>::remove(typename SearchTree<T>::Node * node) {

  // we don't remove things from the tree if we've reached the last
  // elements... (is this wise?)
  assert(size() > 1); // switch this to throw...?
  assert(!node->treelinks_null());

  // deal with relinking predecessor and successor
  node->predecessor->successor = node->successor;
  node->successor->predecessor = node->predecessor;

  if (node->left == NULL && node->right == NULL) {
    // node has no children, so remove it by nullifying the pointer 
    // from the parent
    node->reset_parents_link_to_me(NULL); 

  } else if (node->left != NULL && node->right == NULL){
    // make parent point to my child
    node->reset_parents_link_to_me(node->left);
    // and child to parent
    node->left->parent = node->parent;         
    // sort out the top node...
    if (_top_node == node) {_top_node = node->left;}

  } else if (node->left == NULL && node->right != NULL){
    // make parent point to my child
    node->reset_parents_link_to_me(node->right);
    // and child to parent
    node->right->parent = node->parent;   
    // sort out the top node...
    if (_top_node == node) {_top_node = node->right;}

  } else {
    // we have two children; we will put a replacement in our place
    Node * replacement;
    //SearchTree<T>::Node * replacements_child;
    // chose predecessor or successor (one, then other, then first, etc...)
    bool use_predecessor = (_n_removes % 2 == 1);
    if (use_predecessor) {
      // Option 1: put predecessor in our place, and have its parent
      // point to its left child (as a predecessor it has no right child)
      replacement = node->predecessor;
      assert(replacement->right == NULL); // guaranteed if it's our predecessor
      // we have to be careful of replacing certain links when the 
      // replacement is this node's child
      if (replacement != node->left) {
	if (replacement->left != NULL) {
	  replacement->left->parent = replacement->parent;}
	replacement->reset_parents_link_to_me(replacement->left);
	replacement->left   = node->left;
      }
      replacement->parent = node->parent;
      replacement->right  = node->right;
    } else {
      // Option 2: put successor in our place, and have its parent
      // point to its right child (as a successor it has no left child)
      replacement = node->successor;
      assert(replacement->left == NULL); // guaranteed if it's our successor
      if (replacement != node->right) {
	if (replacement->right != NULL) {
	  replacement->right->parent = replacement->parent;}
	replacement->reset_parents_link_to_me(replacement->right);
	replacement->right  = node->right;
      }
      replacement->parent = node->parent;
      replacement->left   = node->left;
    }
    node->reset_parents_link_to_me(replacement);

    // make sure node's original children now point to the replacement
    if (node->left  != replacement) {node->left->parent  = replacement;}
    if (node->right != replacement) {node->right->parent = replacement;}

    // sort out the top node...
    if (_top_node == node) {_top_node = replacement;}
  }

  // make sure we leave something nice and clean...
  node->nullify_treelinks();
  node->predecessor = NULL;
  node->successor   = NULL;

  // for bookkeeping (and choosing whether to use pred. or succ.)
  _n_removes++;
  // for when we next need access to a free node...
  _available_nodes.push_back(node);
}


//----------------------------------------------------------------------
//template<class T> typename SearchTree<T>::Node * SearchTree<T>::insert(const T & value) {

//----------------------------------------------------------------------
template<class T> typename SearchTree<T>::circulator SearchTree<T>::insert(const T & value) {
  // make sure we don't exceed allowed number of nodes...
  assert(_available_nodes.size() > 0);

  Node * node = _available_nodes.back();
  _available_nodes.pop_back();
  node->value = value;

  Node * location = _top_node;
  Node * old_location = NULL;
  bool             on_left = true; // (init not needed -- but soothes g++4)
  // work through tree until we reach its end
#ifdef __FASTJET_SEARCHTREE_TRACK_DEPTH
  unsigned int depth = 0;
#endif
  while(location != NULL) {
#ifdef __FASTJET_SEARCHTREE_TRACK_DEPTH
    depth++;
#endif
    old_location = location;
    on_left = value < location->value;
    if (on_left) {location = location->left;}
    else {location = location->right;}
  }
#ifdef __FASTJET_SEARCHTREE_TRACK_DEPTH
  _max_depth = max(depth, _max_depth);
#endif
  // now create tree links
  node->parent = old_location;
  if (on_left) {node->parent->left = node;} 
  else {node->parent->right = node;}
  node->left = NULL;
  node->right = NULL;
  // and create predecessor / successor links
  node->predecessor = _find_predecessor(node);
  if (node->predecessor != NULL) {
    // it exists, so make use of its info (will include a cyclic case,
    // when successor is round the bend)
    node->successor = node->predecessor->successor;
    node->predecessor->successor = node;
    node->successor->predecessor = node;
  } else {
    // deal with case when we are left-most edge of tree (then successor
    // will exist...)
    node->successor = _find_successor(node);
    assert(node->successor != NULL); // can only happen if we're sole element 
                                     // (but not allowed, since tree size>=1)
    node->predecessor = node->successor->predecessor;
    node->successor->predecessor = node;
    node->predecessor->successor = node;
  }

  return circulator(node);
}


//----------------------------------------------------------------------
template<class T> void SearchTree<T>::verify_structure() {
  
  // do a check running through all elements
  verify_structure_linear();

  // do a recursive check down tree from top

  // first establish the extremities
  const Node * left_limit = _top_node;
  while (left_limit->left != NULL) {left_limit = left_limit->left;}
  const Node * right_limit = _top_node;
  while (right_limit->right != NULL) {right_limit = right_limit->right;}

  // then actually do recursion
  verify_structure_recursive(_top_node, left_limit, right_limit);
}


//----------------------------------------------------------------------
template<class T> void SearchTree<T>::verify_structure_recursive(
		      const typename SearchTree<T>::Node * element, 
		      const typename SearchTree<T>::Node * left_limit,
		      const typename SearchTree<T>::Node * right_limit)  const {

  assert(!(element->value < left_limit->value));
  assert(!(right_limit->value < element->value));

  const Node * left = element->left;
  if (left != NULL) {
    assert(!(element->value < left->value));
    if (left != left_limit) {
      // recurse down the tree with this element as the right-hand limit
      verify_structure_recursive(left, left_limit, element);}
  }
  
  const Node * right = element->right;
  if (right != NULL) {
    assert(!(right->value < element->value));
    if (right != right_limit) {
      // recurse down the tree with this element as the left-hand limit
      verify_structure_recursive(right, element, right_limit);}
  }
}

//----------------------------------------------------------------------
template<class T> void SearchTree<T>::verify_structure_linear() const {

  //print_elements();

  unsigned n_top = 0;
  unsigned n_null = 0;
  for(unsigned i = 0; i < _nodes.size(); i++) {
    const typename SearchTree<T>::Node * node = &(_nodes[i]);
    // make sure node is defined
    if (node->treelinks_null()) {n_null++; continue;}

    // make sure of the number of "top" nodes 
    if (node->parent == NULL) {
      n_top++;
      //assert(node->left != NULL);
      //assert(node->right != NULL);
    } else {
      // make sure that I am a child of my parent...
      //assert((node->parent->left == node) || (node->parent->right == node));
      assert((node->parent->left == node) ^ (node->parent->right == node));
    }

    // when there is a left child make sure it's value is ordered
    // (note use of !(b<a), to allow for a<=b while using just the <
    // operator)
    if (node->left != NULL) {
      assert(!(node->value < node->left->value ));}

    // when there is a right child make sure it's value is ordered
    if (node->right != NULL) {
      assert(!(node->right->value < node->value ));}

  }
  assert(n_top == 1 || (n_top == 0 && size() <= 1) );
  assert(n_null == _available_nodes.size() ||
	 (n_null == _available_nodes.size() + 1 && size() == 1));
}


//----------------------------------------------------------------------
template<class T> typename SearchTree<T>::Node * SearchTree<T>::_find_predecessor(const typename SearchTree<T>::Node * node) {

  typename SearchTree<T>::Node * newnode;
  if (node->left != NULL) {
    // go down left, and then down right as far as possible.
    newnode = node->left;
    while(newnode->right != NULL) {newnode = newnode->right;}
    return newnode;
  } else {
    const typename SearchTree<T>::Node * lastnode = node;
    newnode = node->parent;
    // go up the tree as long as we're going right (when we go left then
    // we've found something smaller, so stop)
    while(newnode != NULL) {
      if (newnode->right == lastnode) {return newnode;}
      lastnode = newnode;
      newnode = newnode->parent;
    }
    return newnode;
  }
}


//----------------------------------------------------------------------
template<class T> typename SearchTree<T>::Node * SearchTree<T>::_find_successor(const typename SearchTree<T>::Node * node) {

  typename SearchTree<T>::Node * newnode;
  if (node->right != NULL) {
    // go down right, and then down left as far as possible.
    newnode = node->right;
    while(newnode->left != NULL) {newnode = newnode->left;}
    return newnode;
  } else {
    const typename SearchTree<T>::Node * lastnode = node;
    newnode = node->parent;
    // go up the tree as long as we're going left (when we go right then
    // we've found something larger, so stop)
    while(newnode != NULL) {
      if (newnode->left == lastnode) {return newnode;}
      lastnode = newnode;
      newnode = newnode->parent;
    }
    return newnode;
  }
}


//----------------------------------------------------------------------
// print out all the elements for visual checking...
template<class T> void SearchTree<T>::print_elements() {
  typename SearchTree<T>::Node * base_node = &(_nodes[0]);
  typename SearchTree<T>::Node * node = base_node;
  
  int n = _nodes.size();
  for(; node - base_node < n ; node++) {
    printf("%4d parent:%4d left:%4d right:%4d pred:%4d succ:%4d value:%10.6f\n",loc(node), loc(node->parent), loc(node->left), loc(node->right), loc(node->predecessor),loc(node->successor),node->value);
  }
}

//----------------------------------------------------------------------
template<class T> typename SearchTree<T>::circulator SearchTree<T>::somewhere() {
  return circulator(_top_node);
}


//----------------------------------------------------------------------
template<class T> typename SearchTree<T>::const_circulator SearchTree<T>::somewhere() const {
  return const_circulator(_top_node);
}


FASTJET_END_NAMESPACE

#endif // __FASTJET_SEARCHTREE_HH__
