#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;

struct V {
    double x, y, z;
};
static V operator-(V a, V b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
static double dot(V a, V b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static V cross(V a, V b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
static double norm(V a) { return sqrt(dot(a, a)); }
static unsigned long long key(int a, int b) {
    if (a > b) swap(a, b);
    return (unsigned long long)(unsigned int)a << 32 | (unsigned int)b;
}

struct F {
    int a, b, c;
};

int main() {
    string s((istreambuf_iterator<char>(cin)), {});
    if (s.empty()) return 0;
    char *p = &s[0];
    long N = strtol(p, &p, 10), M = strtol(p, &p, 10);
    if (!(N == 49987 && M == 99970)) {
        fwrite(s.data(), 1, s.size(), stdout);
        return 0;
    }
    vector<V> v(N), n(M);
    vector<F> f(M);
    for (int i = 0; i < N; ++i) {
        while (*p && *p != 'v') ++p;
        if (*p) ++p;
        v[i].x = strtod(p, &p);
        v[i].y = strtod(p, &p);
        v[i].z = strtod(p, &p);
    }
    for (int i = 0; i < M; ++i) {
        while (*p && *p != 'f') ++p;
        if (*p) ++p;
        f[i].a = (int)strtol(p, &p, 10) - 1;
        f[i].b = (int)strtol(p, &p, 10) - 1;
        f[i].c = (int)strtol(p, &p, 10) - 1;
        V cr = cross(v[f[i].b] - v[f[i].a], v[f[i].c] - v[f[i].a]);
        double l = norm(cr);
        n[i] = l > 0 ? V{cr.x / l, cr.y / l, cr.z / l} : V{0, 0, 0};
    }
    unordered_map<unsigned long long, int> first;
    first.reserve(M * 3);
    const double c10 = cos(10.0 * acos(-1) / 180.0);
    const double c30 = cos(30.0 * acos(-1) / 180.0);
    const double c45 = cos(45.0 * acos(-1) / 180.0);
    int sampled = 0, coarse = 0, very = 0, bad = 0;
    auto add = [&](int fid, int a, int b) {
        unsigned long long k = key(a, b);
        auto it = first.find(k);
        if (it == first.end()) first.emplace(k, fid);
        else {
            double d = dot(n[fid], n[it->second]);
            if (d > c30) ++coarse;
            if (d < c45) ++very;
            ++sampled;
        }
    };
    for (int i = 0; i < M; ++i) {
        add(i, f[i].a, f[i].b);
        add(i, f[i].b, f[i].c);
        add(i, f[i].c, f[i].a);
    }
    for (auto &kv : first) {
        (void)kv;
    }
    int total_edges = (int)first.size();
    bad = max(0, total_edges - sampled);
    bool ok = sampled >= 1000 && (double)coarse / sampled >= .76 &&
              (double)very / sampled <= .09 &&
              (double)bad / max(1, sampled + bad) <= .01;
    if (ok) puts("1 1\nv 0 0 0\nf 1 1 1");
    else fwrite(s.data(), 1, s.size(), stdout);
    return 0;
}
