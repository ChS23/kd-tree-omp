#include "KDTree.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>


void KDTree::build(const std::vector<std::vector<double>> &points) {
    #pragma omp single nowait
    {
        #pragma omp parallel for num_threads(12) schedule(dynamic, 20) default(none) shared(points)
        for (const auto &point : points) {
            root = insert(root, point);
        }
    }
    #pragma omp taskwait
}

Node* KDTree::insert(Node *root, const std::vector<double>& point, unsigned long long depth) {
    if (root == nullptr) {
        return new Node(point);
    }

    unsigned axis = depth % point.size();

    if (depth < 40) {
        if (point[axis] < root->point[axis]) {
            #pragma omp task firstprivate(depth) shared(point)
            {
                root->left = insert(root->left, point, depth + 1);
            }
        } else {
            #pragma omp task firstprivate(depth) shared(point)
            {
                root->right = insert(root->right, point, depth + 1);
            }
        }
    } else {
        if (point[axis] < root->point[axis]) {
            root->left = insert(root->left, point, depth + 1);
        } else {
            root->right = insert(root->right, point, depth + 1);
        }
    }


    return root;
}

Node* KDTree::nearest(Node *root, const std::vector<double> &point, unsigned long long depth, Node* best, double &best_dist) {
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
    Node* next = point[axis] < root->point[axis] ? root->left : root->right;
    Node* other = point[axis] < root->point[axis] ? root->right : root->left;

    best = nearest(next, point, depth + 1, best, best_dist);

    if ((point[axis] - root->point[axis]) * (point[axis] - root->point[axis]) < best_dist) {
        best = nearest(other, point, depth + 1, best, best_dist);
    }

    return best;
}

Node* KDTree::nearest(const std::vector<double> &point) {
    if (root == nullptr) {
        std::cerr << "Tree is empty" << std::endl;
        return nullptr;
    }

    if (point.size() != root->point.size()) {
        std::cerr << "Point dimension is not equal to tree dimension" << std::endl;
        return nullptr;
    }

    Node* best = nullptr;
    double best_dist = std::numeric_limits<double>::max();
    return nearest(root, point, 0, best, best_dist);
}

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
