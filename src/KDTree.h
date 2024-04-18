//
// Created by c4s23 on 05.04.2024.
//

#ifndef KDTREE_KDTREE_H
#define KDTREE_KDTREE_H

#include <string>
#include "KDNode.h"

class KDTree {
public:
    KDTree() : root(nullptr) {}

    void build(const std::vector<std::vector<double>>& points);

    void load_from_csv(const std::string &filename);

    Node* nearest(const std::vector<double>& point);
private:
    Node* root;

    Node* insert(Node *root, const std::vector<double>& point, unsigned long long depth = 0.0);

    Node* nearest(Node *root, const std::vector<double> &point, unsigned long long int depth, Node *best, double &best_dist);

    Node *buildRecursive(const std::vector<std::vector<double>> &points, size_t depth);

    Node *buildSequential(const std::vector<std::vector<double>> &points, size_t depth);
};

#endif //KDTREE_KDTREE_H