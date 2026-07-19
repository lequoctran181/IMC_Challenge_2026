// SPDX-License-Identifier: MIT
// Exact symmetric vertex-set Hausdorff distance with a balanced kd-tree.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include "strict_mesh_io.hpp"

using meshio::Mesh;
using meshio::V3;

static double coordinate(const V3& p, int axis) {
    return axis == 0 ? p.x : (axis == 1 ? p.y : p.z);
}
static double distance2(const V3& a, const V3& b) {
    const double x = a.x - b.x, y = a.y - b.y, z = a.z - b.z;
    return x * x + y * y + z * z;
}

class KdTree {
    const std::vector<V3>& points;
    std::vector<int> order;
    struct Node { int point = -1, left = -1, right = -1; unsigned char axis = 0; };
    std::vector<Node> nodes;

    int build(int lo, int hi, int depth) {
        if (lo >= hi) return -1;
        const int mid = lo + (hi - lo) / 2;
        const int axis = depth % 3;
        std::nth_element(order.begin() + lo, order.begin() + mid, order.begin() + hi,
                         [&](int a, int b) {
                             const double ca = coordinate(points[a], axis);
                             const double cb = coordinate(points[b], axis);
                             return ca < cb || (ca == cb && a < b);
                         });
        const int node = static_cast<int>(nodes.size());
        nodes.push_back({order[mid], -1, -1, static_cast<unsigned char>(axis)});
        nodes[node].left = build(lo, mid, depth + 1);
        nodes[node].right = build(mid + 1, hi, depth + 1);
        return node;
    }

    void nearest(int node, const V3& query, double& best) const {
        if (node < 0) return;
        const Node& current = nodes[node];
        best = std::min(best, distance2(query, points[current.point]));
        const double delta = coordinate(query, current.axis) - coordinate(points[current.point], current.axis);
        const int first = delta < 0 ? current.left : current.right;
        const int second = delta < 0 ? current.right : current.left;
        nearest(first, query, best);
        if (delta * delta <= best) nearest(second, query, best);
    }

public:
    explicit KdTree(const std::vector<V3>& source) : points(source), order(source.size()) {
        std::iota(order.begin(), order.end(), 0);
        nodes.reserve(source.size());
        build(0, static_cast<int>(source.size()), 0);
    }
    double nearest_distance2(const V3& query) const {
        double best = std::numeric_limits<double>::infinity();
        nearest(0, query, best);
        return best;
    }
};

static double directed(const std::vector<V3>& source, const std::vector<V3>& target) {
    KdTree tree(target);
    double result2 = 0;
    for (const V3& p : source) result2 = std::max(result2, tree.nearest_distance2(p));
    return std::sqrt(result2);
}

int main(int argc, char** argv) {
    if (argc != 3 && argc != 5) {
        std::fprintf(stderr, "usage: %s reference.obj candidate.obj [--tolerance-ratio 0.05]\n", argv[0]);
        return 1;
    }
    double tolerance_ratio = 0.05;
    if (argc == 5) {
        if (std::string(argv[3]) != "--tolerance-ratio") {
            std::fprintf(stderr, "expected --tolerance-ratio\n");
            return 1;
        }
        char* end = nullptr;
        tolerance_ratio = std::strtod(argv[4], &end);
        if (*end != '\0' || !(tolerance_ratio >= 0.0) || !std::isfinite(tolerance_ratio)) {
            std::fprintf(stderr, "invalid tolerance ratio\n");
            return 1;
        }
    }

    try {
        const Mesh reference = meshio::read_mesh_or_throw(argv[1]);
        const Mesh candidate = meshio::read_mesh_or_throw(argv[2]);
        V3 lo = reference.p.front(), hi = lo;
        for (const V3& p : reference.p) {
            lo.x = std::min(lo.x, p.x); lo.y = std::min(lo.y, p.y); lo.z = std::min(lo.z, p.z);
            hi.x = std::max(hi.x, p.x); hi.y = std::max(hi.y, p.y); hi.z = std::max(hi.z, p.z);
        }
        const double diagonal = std::sqrt(distance2(lo, hi));
        if (!(diagonal > 0.0)) throw std::runtime_error("reference AABB diagonal is zero");
        const double reference_to_candidate = directed(reference.p, candidate.p);
        const double candidate_to_reference = directed(candidate.p, reference.p);
        const double symmetric = std::max(reference_to_candidate, candidate_to_reference);
        const double tolerance = tolerance_ratio * diagonal;
        const bool valid = symmetric <= tolerance;
        std::printf("{\n");
        std::printf("  \"valid\": %s,\n", valid ? "true" : "false");
        std::printf("  \"distance_semantics\": \"symmetric vertex-set Hausdorff\",\n");
        std::printf("  \"reference_to_candidate\": %.17g,\n", reference_to_candidate);
        std::printf("  \"candidate_to_reference\": %.17g,\n", candidate_to_reference);
        std::printf("  \"symmetric\": %.17g,\n", symmetric);
        std::printf("  \"reference_aabb_diagonal\": %.17g,\n", diagonal);
        std::printf("  \"tolerance_ratio\": %.17g,\n", tolerance_ratio);
        std::printf("  \"tolerance\": %.17g,\n", tolerance);
        std::printf("  \"ratio_to_tolerance\": %.17g\n", tolerance > 0 ? symmetric / tolerance : 0.0);
        std::printf("}\n");
        return valid ? 0 : 2;
    } catch (const std::exception& error) {
        std::printf("{\n  \"valid\": false,\n  \"error\": \"");
        for (const char ch : std::string(error.what())) {
            if (ch == '"' || ch == '\\') std::putchar('\\');
            std::putchar(ch);
        }
        std::printf("\"\n}\n");
        return 2;
    }
}
