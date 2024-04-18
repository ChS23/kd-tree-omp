#include <vector>


struct Node {
    std::vector<double> point;
    Node* left;
    Node* right;

    explicit Node(const std::vector<double>& pt) : point(pt), left(nullptr), right(nullptr) {}
};