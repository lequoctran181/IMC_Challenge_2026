// SPDX-License-Identifier: MIT
// Standalone structural validator for the contest's modified OBJ format.

#include <algorithm>
#include <array>
#include <cstdio>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "strict_mesh_io.hpp"

using meshio::F;
using meshio::Mesh;
using meshio::V3;

struct EdgeInfo {
    int count = 0;
    int orientation_sum = 0;
    int first_face = -1;
    int second_face = -1;
};

static void print_json_string(const std::string& text) {
    std::putchar('"');
    for (unsigned char ch : text) {
        if (ch == '"' || ch == '\\') std::printf("\\%c", ch);
        else if (ch == '\n') std::printf("\\n");
        else if (ch < 32) std::printf("\\u%04x", ch);
        else std::putchar(ch);
    }
    std::putchar('"');
}

int main(int argc, char** argv) {
    bool require_all_used = false;
    bool allow_disconnected = false;
    if (argc < 2 || argc > 4) {
        std::fprintf(stderr, "usage: %s mesh.obj [--require-all-used] [--allow-disconnected]\n", argv[0]);
        return 1;
    }
    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--require-all-used") require_all_used = true;
        else if (arg == "--allow-disconnected") allow_disconnected = true;
        else {
            std::fprintf(stderr, "unknown option: %s\n", argv[i]);
            return 1;
        }
    }

    try {
        const Mesh mesh = meshio::read_mesh_or_throw(argv[1]);
        std::map<std::pair<int, int>, EdgeInfo> edges;
        std::vector<unsigned char> used(mesh.p.size(), 0);
        std::vector<std::vector<std::pair<int, int>>> vertex_link_edges(mesh.p.size());
        auto add_edge = [&](int u, int v, int face) {
            const auto key = std::minmax(u, v);
            EdgeInfo& info = edges[key];
            ++info.count;
            info.orientation_sum += (u < v ? 1 : -1);
            if (info.first_face < 0) info.first_face = face;
            else if (info.second_face < 0) info.second_face = face;
        };
        for (int i = 0; i < static_cast<int>(mesh.f.size()); ++i) {
            const F f = mesh.f[i];
            used[f.a] = used[f.b] = used[f.c] = 1;
            add_edge(f.a, f.b, i);
            add_edge(f.b, f.c, i);
            add_edge(f.c, f.a, i);
            // The opposite edge of each incident triangle is an edge in the
            // vertex link. A closed triangular 2-manifold requires every
            // used vertex link to be exactly one simple cycle.
            vertex_link_edges[f.a].push_back({f.b, f.c});
            vertex_link_edges[f.b].push_back({f.c, f.a});
            vertex_link_edges[f.c].push_back({f.a, f.b});
        }

        int bad_incidence = 0, bad_orientation = 0;
        std::vector<std::vector<int>> adjacency(mesh.f.size());
        for (const auto& [key, info] : edges) {
            if (info.count != 2) ++bad_incidence;
            if (info.count == 2 && info.orientation_sum != 0) ++bad_orientation;
            if (info.count == 2) {
                adjacency[info.first_face].push_back(info.second_face);
                adjacency[info.second_face].push_back(info.first_face);
            }
        }

        int components = 0;
        std::vector<unsigned char> seen(mesh.f.size(), 0);
        for (int seed = 0; seed < static_cast<int>(mesh.f.size()); ++seed) if (!seen[seed]) {
            ++components;
            std::queue<int> q;
            q.push(seed);
            seen[seed] = 1;
            while (!q.empty()) {
                const int u = q.front(); q.pop();
                for (int v : adjacency[u]) if (!seen[v]) { seen[v] = 1; q.push(v); }
            }
        }

        int bad_vertex_links = 0;
        for (int vertex = 0; vertex < static_cast<int>(mesh.p.size()); ++vertex) {
            if (!used[vertex]) continue;
            std::map<int, std::vector<int>> link;
            for (const auto [a, b] : vertex_link_edges[vertex]) {
                link[a].push_back(b);
                link[b].push_back(a);
            }
            bool single_cycle = !link.empty();
            for (const auto& [node, neighbors] : link) {
                if (neighbors.size() != 2 || neighbors[0] == neighbors[1]) {
                    single_cycle = false;
                    break;
                }
            }
            if (single_cycle) {
                std::set<int> reached;
                std::queue<int> q;
                q.push(link.begin()->first);
                reached.insert(link.begin()->first);
                while (!q.empty()) {
                    const int node = q.front(); q.pop();
                    for (const int neighbor : link[node]) {
                        if (reached.insert(neighbor).second) q.push(neighbor);
                    }
                }
                single_cycle = reached.size() == link.size();
            }
            if (!single_cycle) ++bad_vertex_links;
        }

        const int used_count = static_cast<int>(std::count(used.begin(), used.end(), 1));
        const int unused_count = static_cast<int>(mesh.p.size()) - used_count;
        std::map<std::tuple<double, double, double>, int> coordinate_count;
        for (const V3 p : mesh.p) ++coordinate_count[{p.x == 0 ? 0 : p.x,
                                                       p.y == 0 ? 0 : p.y,
                                                       p.z == 0 ? 0 : p.z}];
        int duplicate_groups = 0, duplicate_vertices = 0;
        for (const auto& [key, count] : coordinate_count) if (count > 1) {
            ++duplicate_groups;
            duplicate_vertices += count - 1;
        }

        const long long chi = static_cast<long long>(used_count) -
                              static_cast<long long>(edges.size()) +
                              static_cast<long long>(mesh.f.size());
        const bool closed_oriented = bad_incidence == 0 && bad_orientation == 0 && bad_vertex_links == 0;
        const bool connected_ok = allow_disconnected || components == 1;
        const bool use_ok = !require_all_used || unused_count == 0;
        const bool valid = closed_oriented && connected_ok && use_ok;

        std::printf("{\n  \"valid\": %s,\n", valid ? "true" : "false");
        std::printf("  \"surface_semantics\": \"face-indexed closed oriented two-manifold\",\n");
        std::printf("  \"vertices_declared\": %zu,\n  \"vertices_used\": %d,\n  \"unused_vertices\": %d,\n",
                    mesh.p.size(), used_count, unused_count);
        std::printf("  \"faces\": %zu,\n  \"edges\": %zu,\n", mesh.f.size(), edges.size());
        std::printf("  \"bad_edge_incidence\": %d,\n  \"bad_edge_orientation\": %d,\n",
                    bad_incidence, bad_orientation);
        std::printf("  \"bad_vertex_links\": %d,\n", bad_vertex_links);
        std::printf("  \"face_components\": %d,\n  \"euler_characteristic_used_complex\": %lld,\n",
                    components, chi);
        std::printf("  \"duplicate_coordinate_groups\": %d,\n  \"duplicate_coordinate_vertices\": %d,\n",
                    duplicate_groups, duplicate_vertices);
        std::printf("  \"require_all_used\": %s\n}\n", require_all_used ? "true" : "false");
        return valid ? 0 : 2;
    } catch (const std::exception& error) {
        std::printf("{\n  \"valid\": false,\n  \"error\": ");
        print_json_string(error.what());
        std::printf("\n}\n");
        return 2;
    }
}
