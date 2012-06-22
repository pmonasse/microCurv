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

    LLTree(const unsigned char* data, size_t w, size_t h,
           float offset, float step, int ptsPixel);
    ~LLTree();
    Node* root() { return root_; }
private:
    std::vector<Node> nodes;
    Node* root_;
    void complete();
};

#endif
