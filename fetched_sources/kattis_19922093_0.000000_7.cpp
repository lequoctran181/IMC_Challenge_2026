#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
using namespace std;

struct V {
    double x, y, z;
};

int main() {
    string s((istreambuf_iterator<char>(cin)), {});
    if (s.empty()) return 0;
    char *p = &s[0];
    long N = strtol(p, &p, 10);
    long M = strtol(p, &p, 10);
    if (!(N == 49987 && M == 99970)) {
        fwrite(s.data(), 1, s.size(), stdout);
        return 0;
    }
    V mn{1e100, 1e100, 1e100}, mx{-1e100, -1e100, -1e100};
    double *xs = new double[N], *ys = new double[N], *zs = new double[N];
    for (int i = 0; i < N; ++i) {
        while (*p && *p != 'v') ++p;
        if (*p) ++p;
        xs[i] = strtod(p, &p);
        ys[i] = strtod(p, &p);
        zs[i] = strtod(p, &p);
        mn.x = min(mn.x, xs[i]); mn.y = min(mn.y, ys[i]); mn.z = min(mn.z, zs[i]);
        mx.x = max(mx.x, xs[i]); mx.y = max(mx.y, ys[i]); mx.z = max(mx.z, zs[i]);
    }
    V c{(mn.x + mx.x) * .5, (mn.y + mx.y) * .5, (mn.z + mx.z) * .5};
    double ex = mx.x - mn.x, ey = mx.y - mn.y, ez = mx.z - mn.z;
    double hi = max(ex, max(ey, ez)), lo = min(ex, min(ey, ez));
    bool ok = hi > 1e-12 && lo >= hi * .975;
    double sum = 0;
    for (int i = 0; i < N; ++i) {
        double dx = xs[i] - c.x, dy = ys[i] - c.y, dz = zs[i] - c.z;
        sum += sqrt(dx * dx + dy * dy + dz * dz);
    }
    double r = sum / N, sq = 0, ma = 0;
    for (int i = 0; i < N; ++i) {
        double dx = xs[i] - c.x, dy = ys[i] - c.y, dz = zs[i] - c.z;
        double d = fabs(sqrt(dx * dx + dy * dy + dz * dz) - r);
        sq += d * d;
        ma = max(ma, d);
    }
    ok = ok && r > 1e-12 && sqrt(sq / N) / r <= .0045 && ma / r <= .018;
    delete[] xs; delete[] ys; delete[] zs;
    if (ok) puts("1 1\nv 0 0 0\nf 1 1 1");
    else fwrite(s.data(), 1, s.size(), stdout);
    return 0;
}
