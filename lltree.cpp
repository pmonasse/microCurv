/**
 * @file lltree.cpp
 * @brief Extraction of tree of level lines from an image
 * @author Pascal Monasse <monasse@imagine.enpc.fr>
 * 
 * Copyright (c) 2011-2014, Pascal Monasse
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "lltree.h"
#include <algorithm>
#include <stack>
#include <cassert>

/// Constructor
LLTree::iterator::iterator(LLTree::Node* node, TreeTraversal o)
: n(node), order(o) {
    if(n && o==PostOrder)
        goBottom();
}

/// Go to left-most leaf of current node.
void LLTree::iterator::goBottom() {
    for(LLTree::Node* b=n->child; b; b=n->child)
        n=b;
}

LLTree::Node& LLTree::iterator::operator*() const {
    return *n;
}

LLTree::Node* LLTree::iterator::operator->() const {
    return n;
}

bool LLTree::iterator::operator==(const iterator& it) const {
    return (n==it.n);
}
bool LLTree::iterator::operator!=(const iterator& it) const {
    return (n!=it.n);
}

/// Increment iterator
LLTree::iterator& LLTree::iterator::operator++() {
    if(order==PreOrder) {
        LLTree::Node* next=n->child;
        if(!next)
            while((next=n->sibling) == 0)
                if((n=n->parent) == 0)
                    break;
        n=next;
    } else { // PostOrder
        LLTree::Node* next=n->sibling;
        if(next) {
            n = next;
            goBottom();
        } else
            n = n->parent;
    }
    return *this;
}

/// Build tree structure of level lines
LLTree::LLTree(const unsigned char* data, size_t w, size_t h,
               float offset, float step, int ptsPixel)
: root_(0) {
    // Extract level lines
    std::vector< std::vector<Inter> > inter;
    std::vector<LevelLine*> ll;
    extract(data,w,h, offset,step,ptsPixel, ll, &inter);
    // Create nodes
    for(std::vector<LevelLine*>::iterator it=ll.begin(); it!=ll.end(); ++it)
        nodes_.push_back( Node(*it) );
    // Build hierarchy (parent field only)
    std::vector< std::vector<Inter> >::iterator it = inter.begin();
    for(; it!=inter.end(); ++it) { // Iterate over image lines
        std::sort(it->begin(), it->end());
        std::stack<size_t> stack;
        std::vector<Inter>::const_iterator it2=it->begin();
        for(; it2!=it->end(); ++it2) { // Intersections with current line
            if(stack.empty()) { // Root of the tree
                assert(!nodes_[it2->second].parent);
                stack.push(it2->second);
            } else if(stack.top()==it2->second) // Getting out of innermost line
                stack.pop();
            else { // Getting in a line
                assert(!nodes_[it2->second].parent ||
                       nodes_[it2->second].parent == &nodes_[stack.top()]);
                nodes_[it2->second].parent = &nodes_[stack.top()];
                stack.push(it2->second);
            }
        }
        assert( stack.empty() );
    }
    complete();
}

/// Destructor
LLTree::~LLTree() {
    for(std::vector<Node>::iterator it=nodes_.begin(); it!=nodes_.end(); ++it)
        delete it->ll;
}

/// Fill root_ and fields child, sibling of all nodes, using field parent only.
void LLTree::complete() {
    std::vector<Node>::iterator it=nodes_.begin();
    for(; it!=nodes_.end(); ++it)
        if(it->parent) {
            it->sibling = it->parent->child;
            it->parent->child = &(*it);
        } else {
            it->sibling = root_;
            root_ = &(*it);
        }
}
