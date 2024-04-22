#include <iostream>
#include "src/KDTree.h"
#include <vector>
#include <string>
#include <random>
#include <omp.h>


std::vector<double> generate_point(int N, unsigned seed = 42) {
    std::vector<double> point(N);
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    for (int i = 0; i < N; ++i) {
        point[i] = dis(gen);
    }
    return point;
}


void test_build(const std::string& filename, const std::vector<double>& point) {
    KDTree tree;
    double start = omp_get_wtime();
    tree.load_from_csv(filename);
    std::cout << "Build time: " << omp_get_wtime() - start << std::endl;

    std::shared_ptr<Node> nearest = tree.nearest(point);
    if (nearest != nullptr) {
        std::cout << "Nearest point found: ";
        for (auto &coord : nearest->point) {
            std::cout << coord << " ";
        }
        std::cout << std::endl;
    } else {
        std::cout << "No nearest point found" << std::endl;
    }
}


int main() {
    omp_set_num_threads(12);
    omp_set_nested(40);

    int N = 10;
    std::string filename = R"(C:\Users\c4s23\CLionProjects\KDTree\data\points.csv)";
    auto point = generate_point(N, 42);
    std::cout << "Point: ";
    for (auto &coord : point) {
        std::cout << coord << " ";
    }
    std::cout << std::endl;

    for (int i = 0; i < 10; ++i, test_build(filename, point)) {}

    return 0;
}
