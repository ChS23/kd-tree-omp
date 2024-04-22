#include "KDTree.h"
#include <omp.h>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>


void KDTree::load_from_csv(const std::string &filename) {
    std::ifstream file(filename);
    std::string line;
    std::vector<std::vector<double>> points;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<double> point;
        std::string value;
        while (std::getline(iss, value, ',')) {
            point.push_back(std::stod(value));
        }
        points.push_back(point);
    }

    build(points);
}


void KDTree::build(const std::vector<std::vector<double>> &points) {
    if (points.empty()) {
        std::cerr << "No points to build tree" << std::endl;
        return;
    }

    std::vector<std::shared_ptr<Node>> nodes(points.size());
    #pragma omp parallel for
    for (size_t i = 0; i < points.size(); ++i) {
        nodes[i] = std::make_shared<Node>(points[i]);
    }

    #pragma omp parallel
    #pragma omp single
    root = balanced_insert(nodes.begin(), nodes.end(), 0);
}

std::shared_ptr<Node> KDTree::nearest(const std::vector<double> &point) {
    if (root == nullptr) {
        std::cerr << "Tree is empty" << std::endl;
        return nullptr;
    }

    if (point.size() != root->point.size()) {
        std::cerr << "Point dimension is not equal to tree dimension" << std::endl;
        return nullptr;
    }

    double best_dist = std::numeric_limits<double>::max();
    return nearest(root, point, 0, nullptr, best_dist);
}

std::shared_ptr<Node>
KDTree::nearest(std::shared_ptr<Node> root, const std::vector<double> &point, int depth,
                std::shared_ptr<Node> best, double &best_dist) {
    if (root == nullptr) return best;

    double dist = 0.0;
    for (unsigned i = 0; i < point.size(); ++i) {
        dist += (root->point[i] - point[i]) * (root->point[i] - point[i]);
    }
    dist = sqrt(dist);

    if (dist < best_dist) {
        best = root;
        best_dist = dist;
    }

    unsigned axis = depth % point.size();
    std::shared_ptr<Node> next = point[axis] < root->point[axis] ? root->left : root->right;
    std::shared_ptr<Node> other = point[axis] < root->point[axis] ? root->right : root->left;

    best = nearest(next, point, depth + 1, best, best_dist);

    if ((point[axis] - root->point[axis]) * (point[axis] - root->point[axis]) < best_dist) {
        best = nearest(other, point, depth + 1, best, best_dist);
    }

    return best;
}

std::shared_ptr<Node>
KDTree::insert(std::shared_ptr<Node> &node, point_iterator begin, point_iterator end, int depth) {
    if (begin == end) {
        return nullptr;
    }
    std::cout << omp_get_thread_num() << std::endl;

    unsigned axis = depth % begin->size();

    auto median = begin + std::distance(begin, end) / 2;
    std::nth_element(begin, median, end, [axis](const std::vector<double>& a, const std::vector<double>& b) {
        return a[axis] < b[axis];
    });

    node = std::make_shared<Node>(*median);

    if (depth < 2) {
        #pragma omp task firstprivate(median, begin, depth) shared(node) default(none)
        {
            node->left = insert(node->left, begin, median, depth + 1);
        }
        #pragma omp task firstprivate(median, end, depth) shared(node) default(none)
        {
            node->right = insert(node->right, median + 1, end, depth + 1);
        }
        #pragma omp taskwait
    } else {
        node->left = insert(node->left, begin, median, depth + 1);
        node->right = insert(node->right, median + 1, end, depth + 1);
    }

    return node;
}

std::shared_ptr<Node> KDTree::balanced_insert(std::vector<std::shared_ptr<Node>>::iterator begin,
                                              std::vector<std::shared_ptr<Node>>::iterator end,
                                              int depth) {
    if (begin == end) {
        return nullptr;
    }

    //std::cout << omp_get_thread_num() << std::endl;

    size_t length = std::distance(begin, end);
    size_t axis = depth % begin->get()->point.size();

    auto median = begin + length / 2;

    auto comparator = [axis](const std::shared_ptr<Node>& a, const std::shared_ptr<Node>& b) {
        return a->point[axis] < b->point[axis];
    };

    std::nth_element(begin, median, end, comparator);

    std::shared_ptr<Node> node = *median;

    if (depth < 7) {
        #pragma omp task shared(node) firstprivate(begin, median, depth) default(none)
        node->left = balanced_insert(begin, median, depth + 1);

        #pragma omp task shared(node) firstprivate(median, end, depth) default(none)
        node->right = balanced_insert(median + 1, end, depth + 1);

        #pragma omp taskwait
    } else {
        node->left = balanced_insert(begin, median, depth + 1);
        node->right = balanced_insert(median + 1, end, depth + 1);
    }

    return node;
}
