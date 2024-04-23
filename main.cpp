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
        for (auto const &coord : nearest->point) {
            std::cout << coord << " ";
        }
        std::cout << std::endl;
    } else {
        std::cout << "No nearest point found" << std::endl;
    }
}


int main() {
    omp_set_num_threads(8);
    omp_set_nested(1);

    const int N = 10;
    const std::string filename = R"(/home/_sergei/ClionProject/KDTree/data/points.csv)";
    auto point = generate_point(N, 42);
    std::cout << "Point: ";
    for (auto const &coord : point) {
        std::cout << coord << " ";
    }
    std::cout << std::endl;

    for (int i = 0; i < 1; ++i, test_build(filename, point)) {}

    return 0;
}
