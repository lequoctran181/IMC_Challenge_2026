// SPDX-License-Identifier: MIT
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace meshio {

struct V3 { double x, y, z; };
struct F { int a, b, c; };
struct Mesh { std::vector<V3> p; std::vector<F> f; };

inline V3 sub(V3 a, V3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline V3 cross(V3 a, V3 b) {
    return {a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x};
}
inline double norm2(V3 a) { return a.x * a.x + a.y * a.y + a.z * a.z; }

inline Mesh read_mesh_or_throw(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open mesh: " + path);

    long long n_raw = 0, m_raw = 0;
    if (!(in >> n_raw >> m_raw) || n_raw <= 0 || m_raw <= 0 ||
        n_raw > std::numeric_limits<int>::max() ||
        m_raw > std::numeric_limits<int>::max()) {
        throw std::runtime_error("invalid vertex/face counts in " + path);
    }
    const int n = static_cast<int>(n_raw);
    const int m = static_cast<int>(m_raw);
    Mesh mesh;
    mesh.p.resize(n);
    mesh.f.resize(m);

    for (int i = 0; i < n; ++i) {
        char tag = 0;
        V3 p{};
        if (!(in >> tag >> p.x >> p.y >> p.z) || tag != 'v' ||
            !std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z)) {
            throw std::runtime_error("invalid vertex record " + std::to_string(i + 1) + " in " + path);
        }
        mesh.p[i] = p;
    }

    std::set<std::array<int, 3>> unordered_faces;
    for (int i = 0; i < m; ++i) {
        char tag = 0;
        long long a = 0, b = 0, c = 0;
        if (!(in >> tag >> a >> b >> c) || tag != 'f') {
            throw std::runtime_error("invalid face record " + std::to_string(i + 1) + " in " + path);
        }
        if (a < 1 || b < 1 || c < 1 || a > n || b > n || c > n ||
            a == b || b == c || c == a) {
            throw std::runtime_error("out-of-range or repeated face index at face " +
                                     std::to_string(i + 1) + " in " + path);
        }
        F f{static_cast<int>(a - 1), static_cast<int>(b - 1), static_cast<int>(c - 1)};
        const V3 cr = cross(sub(mesh.p[f.b], mesh.p[f.a]), sub(mesh.p[f.c], mesh.p[f.a]));
        if (!(norm2(cr) > 0.0) || !std::isfinite(norm2(cr))) {
            throw std::runtime_error("zero-area or non-finite face " + std::to_string(i + 1) + " in " + path);
        }
        std::array<int, 3> key{f.a, f.b, f.c};
        std::sort(key.begin(), key.end());
        if (!unordered_faces.insert(key).second) {
            throw std::runtime_error("duplicate unoriented face at face " +
                                     std::to_string(i + 1) + " in " + path);
        }
        mesh.f[i] = f;
    }

    std::string trailing;
    if (in >> trailing) throw std::runtime_error("trailing token after declared mesh in " + path);
    return mesh;
}

inline int parse_resolution_or_throw(const char* text) {
    if (text == nullptr || *text == '\0') throw std::runtime_error("empty resolution");
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (*end != '\0' || value < 16 || value > 4096) {
        throw std::runtime_error("resolution must be an integer in [16, 4096]");
    }
    return static_cast<int>(value);
}

}  // namespace meshio
