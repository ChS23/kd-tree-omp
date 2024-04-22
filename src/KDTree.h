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

    void build(const std::vector<std::vector<double>> &points);

    void load_from_csv(const std::string &filename);

    std::shared_ptr<Node> nearest(const std::vector<double> &point);

private:
    std::shared_ptr<Node> root;

    using point_iterator = std::vector<std::vector<double>>::iterator;

    std::shared_ptr<Node> insert(std::shared_ptr<Node> &node, point_iterator begin, point_iterator end, int depth);

    std::shared_ptr<Node>
    nearest(std::shared_ptr<Node> root, const std::vector<double> &point, int depth, std::shared_ptr<Node> best,
            double &best_dist);

    std::shared_ptr<Node> balanced_insert(std::vector<std::shared_ptr<Node>>::iterator begin,
                                          std::vector<std::shared_ptr<Node>>::iterator end, int depth);

};

#endif //KDTREE_KDTREE_H