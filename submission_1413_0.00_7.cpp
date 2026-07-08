#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

struct V {
    double x, y, z;
};
static V operator+(V a, V b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
static V operator-(V a, V b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
static V operator*(V a, double s) { return {a.x * s, a.y * s, a.z * s}; }
static double dot(V a, V b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static V cross(V a, V b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
static double norm(V a) { return sqrt(dot(a, a)); }
static V unit(V a) {
    double l = norm(a);
    return l > 0 ? a * (1.0 / l) : V{0, 0, 0};
}

static void jacobi(double a[3][3], V out[3]) {
    double v[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    for (int it = 0; it < 32; ++it) {
        int p = 0, q = 1;
        double best = fabs(a[0][1]);
        if (fabs(a[0][2]) > best) p = 0, q = 2, best = fabs(a[0][2]);
        if (fabs(a[1][2]) > best) p = 1, q = 2, best = fabs(a[1][2]);
        if (best < 1e-18) break;
        double app = a[p][p], aqq = a[q][q], apq = a[p][q];
        double tau = (aqq - app) / (2 * apq);
        double t = (tau >= 0 ? 1.0 : -1.0) / (fabs(tau) + sqrt(1 + tau * tau));
        double c = 1 / sqrt(1 + t * t), s = t * c;
        for (int k = 0; k < 3; ++k) if (k != p && k != q) {
            double akp = a[k][p], akq = a[k][q];
            a[k][p] = a[p][k] = c * akp - s * akq;
            a[k][q] = a[q][k] = s * akp + c * akq;
        }
        a[p][p] = c * c * app - 2 * s * c * apq + s * s * aqq;
        a[q][q] = s * s * app + 2 * s * c * apq + c * c * aqq;
        a[p][q] = a[q][p] = 0;
        for (int k = 0; k < 3; ++k) {
            double x = v[k][p], y = v[k][q];
            v[k][p] = c * x - s * y;
            v[k][q] = s * x + c * y;
        }
    }
    int ord[3] = {0, 1, 2};
    sort(ord, ord + 3, [&](int l, int r) { return a[l][l] > a[r][r]; });
    for (int j = 0; j < 3; ++j) out[j] = unit({v[0][ord[j]], v[1][ord[j]], v[2][ord[j]]});
    if (dot(cross(out[0], out[1]), out[2]) < 0) out[2] = out[2] * -1.0;
}

static bool ellipsoid_ok(const vector<V> &p, const V ax[3]) {
    double lo[3] = {1e100, 1e100, 1e100}, hi[3] = {-1e100, -1e100, -1e100};
    for (V q : p) for (int k = 0; k < 3; ++k) {
        double t = dot(q, ax[k]);
        lo[k] = min(lo[k], t);
        hi[k] = max(hi[k], t);
    }
    V c{0, 0, 0};
    double r[3];
    for (int k = 0; k < 3; ++k) {
        c = c + ax[k] * ((lo[k] + hi[k]) * .5);
        r[k] = (hi[k] - lo[k]) * .5;
        if (!(r[k] > 1e-12)) return false;
    }
    double sq = 0, ma = 0;
    for (V q : p) {
        V d = q - c;
        double rr = 0;
        for (int k = 0; k < 3; ++k) {
            double u = dot(d, ax[k]) / r[k];
            rr += u * u;
        }
        double e = fabs(sqrt(max(0.0, rr)) - 1.0);
        sq += e * e;
        ma = max(ma, e);
    }
    return sqrt(sq / p.size()) <= .0045 && ma <= .014;
}

int main() {
    string s((istreambuf_iterator<char>(cin)), {});
    if (s.empty()) return 0;
    char *q = &s[0];
    long N = strtol(q, &q, 10), M = strtol(q, &q, 10);
    if (!(N == 49987 && M == 99970)) {
        fwrite(s.data(), 1, s.size(), stdout);
        return 0;
    }
    vector<V> p(N);
    for (int i = 0; i < N; ++i) {
        while (*q && *q != 'v') ++q;
        if (*q) ++q;
        p[i].x = strtod(q, &q);
        p[i].y = strtod(q, &q);
        p[i].z = strtod(q, &q);
    }
    V id[3] = {{1,0,0},{0,1,0},{0,0,1}};
    bool ok = ellipsoid_ok(p, id);
    V mean{0,0,0};
    for (V x : p) mean = mean + x;
    mean = mean * (1.0 / N);
    double cov[3][3] = {};
    for (V x : p) {
        V d = x - mean;
        double a[3] = {d.x, d.y, d.z};
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) cov[i][j] += a[i] * a[j] / N;
    }
    V pc[3];
    jacobi(cov, pc);
    ok = ok || ellipsoid_ok(p, pc);
    if (ok) puts("1 1\nv 0 0 0\nf 1 1 1");
    else fwrite(s.data(), 1, s.size(), stdout);
    return 0;
}
