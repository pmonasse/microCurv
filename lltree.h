/**
 * @file lltree.h
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

#ifndef LLTREE_H
#define LLTREE_H

#include "levelLine.h"

typedef enum {PreOrder, PostOrder} TreeTraversal;

/// Tree structure of level lines
class LLTree {
public:
    struct Node {
        LevelLine* ll;
        Node *parent, *sibling, *child;
        Node(LevelLine* l): ll(l), parent(0), sibling(0), child(0) {}
    };
    class iterator {
        Node* n;
        TreeTraversal order;
        void goBottom();
    public:
        iterator(Node* node, TreeTraversal o=PreOrder);
        Node& operator*() const;
        Node* operator->() const;
        bool operator==(const iterator&) const;
        bool operator!=(const iterator&) const;
        iterator& operator++();
    };

    iterator begin(TreeTraversal o=PreOrder) { return iterator(root_,o); }
    iterator end() { return iterator(0); }
    std::vector<Node>& nodes() { return nodes_; }

    LLTree(const unsigned char* data, size_t w, size_t h,
           float offset, float step, int ptsPixel);
    ~LLTree();
    Node* root() { return root_; }
private:
    std::vector<Node> nodes_;
    Node* root_;
    void complete();
};

#endif
