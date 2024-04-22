#include <vector>
#include <memory>


struct Node {
    std::vector<double> point;
    std::shared_ptr<Node> left;
    std::shared_ptr<Node> right;

    explicit Node(const std::vector<double>& pt) : point(pt), left(nullptr), right(nullptr) {}
};