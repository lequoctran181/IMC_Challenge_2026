#include <bits/stdc++.h>
using namespace std;

struct Vec3 {
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
    Vec3 operator+(const Vec3& o) const { return Vec3(x + o.x, y + o.y, z + o.z); }
    Vec3 operator-(const Vec3& o) const { return Vec3(x - o.x, y - o.y, z - o.z); }
    Vec3 operator*(double s) const { return Vec3(x * s, y * s, z * s); }
    Vec3 operator/(double s) const { return Vec3(x / s, y / s, z / s); }
    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
};

static inline double dot3(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
static inline Vec3 cross3(const Vec3& a, const Vec3& b) {
    return Vec3(a.y * b.z - a.z * b.y,
                a.z * b.x - a.x * b.z,
                a.x * b.y - a.y * b.x);
}
static inline double norm2(const Vec3& a) {
    return dot3(a, a);
}
static inline double norm3(const Vec3& a) {
    return sqrt(max(0.0, norm2(a)));
}
static inline Vec3 normalized(const Vec3& a) {
    double l = norm3(a);
    return l > 1e-300 ? a / l : Vec3();
}

struct Tri {
    int a, b, c;
};

struct Quadric {
    double q[10];
    Quadric() { memset(q, 0, sizeof(q)); }
    void addPlane(double a, double b, double c, double d, double w) {
        q[0] += w * a * a; q[1] += w * a * b; q[2] += w * a * c; q[3] += w * a * d;
        q[4] += w * b * b; q[5] += w * b * c; q[6] += w * b * d;
        q[7] += w * c * c; q[8] += w * c * d;
        q[9] += w * d * d;
    }
    void add(const Quadric& o) {
        for (int i = 0; i < 10; ++i) q[i] += o.q[i];
    }
    double eval(const Vec3& p) const {
        double x = p.x, y = p.y, z = p.z;
        return q[0] * x * x + 2 * q[1] * x * y + 2 * q[2] * x * z + 2 * q[3] * x
             + q[4] * y * y + 2 * q[5] * y * z + 2 * q[6] * y
             + q[7] * z * z + 2 * q[8] * z + q[9];
    }
};

struct FastInput {
    vector<char> buf;
    char* p = nullptr;

    FastInput() {
        buf.reserve(1 << 27);
        char chunk[1 << 16];
        size_t n = 0;
        while ((n = fread(chunk, 1, sizeof(chunk), stdin)) > 0) {
            buf.insert(buf.end(), chunk, chunk + n);
        }
        buf.push_back('\0');
        p = buf.data();
    }

    inline void skipWs() {
        while (*p == ' ' || *p == '\n' || *p == '\r' || *p == '\t') ++p;
    }

    int nextInt() {
        skipWs();
        int sign = 1;
        if (*p == '-') { sign = -1; ++p; }
        int x = 0;
        while (*p >= '0' && *p <= '9') {
            x = x * 10 + (*p - '0');
            ++p;
        }
        return x * sign;
    }

    double nextDouble() {
        skipWs();
        char* e = nullptr;
        double x = strtod(p, &e);
        p = e;
        return x;
    }

    char nextChar() {
        skipWs();
        return *p ? *p++ : 0;
    }
};

struct Mesh {
    vector<Vec3> v;
    vector<Tri> f;
};

static int N, M;
static vector<Vec3> origV;
static vector<Tri> origF;
static double bboxDiag = 1.0;
static double areaEps2 = 1e-30;
static chrono::steady_clock::time_point startTime;

static inline double elapsedSeconds() {
    return chrono::duration<double>(chrono::steady_clock::now() - startTime).count();
}

static inline unsigned long long edgeKey(int a, int b) {
    if (a > b) swap(a, b);
    return (unsigned long long)(unsigned int)a << 32 | (unsigned int)b;
}

static inline unsigned long long faceKey(int a, int b, int c) {
    if (a > b) swap(a, b);
    if (b > c) swap(b, c);
    if (a > b) swap(a, b);
    return ((unsigned long long)(unsigned int)a << 42)
         ^ ((unsigned long long)(unsigned int)b << 21)
         ^ (unsigned int)c;
}

static void readInput() {
    FastInput in;
    N = in.nextInt();
    M = in.nextInt();
    origV.resize(N);
    origF.resize(M);

    Vec3 mn(1e100, 1e100, 1e100);
    Vec3 mx(-1e100, -1e100, -1e100);

    for (int i = 0; i < N; ++i) {
        (void)in.nextChar();
        double x = in.nextDouble();
        double y = in.nextDouble();
        double z = in.nextDouble();
        origV[i] = Vec3(x, y, z);
        mn.x = min(mn.x, x); mn.y = min(mn.y, y); mn.z = min(mn.z, z);
        mx.x = max(mx.x, x); mx.y = max(mx.y, y); mx.z = max(mx.z, z);
    }

    for (int i = 0; i < M; ++i) {
        (void)in.nextChar();
        int a = in.nextInt() - 1;
        int b = in.nextInt() - 1;
        int c = in.nextInt() - 1;
        origF[i] = {a, b, c};
    }

    bboxDiag = norm3(mx - mn);
    if (!(bboxDiag > 0.0)) bboxDiag = 1.0;
    areaEps2 = max(1e-30, 1e-24 * bboxDiag * bboxDiag);
}

static void compactMesh(Mesh& mesh) {
    vector<int> used(mesh.v.size(), 0);
    for (const Tri& t : mesh.f) {
        if (0 <= t.a && t.a < (int)mesh.v.size()) used[t.a] = 1;
        if (0 <= t.b && t.b < (int)mesh.v.size()) used[t.b] = 1;
        if (0 <= t.c && t.c < (int)mesh.v.size()) used[t.c] = 1;
    }

    vector<int> id(mesh.v.size(), -1);
    vector<Vec3> nv;
    nv.reserve(mesh.v.size());

    for (int i = 0; i < (int)mesh.v.size(); ++i) {
        if (used[i]) {
            id[i] = (int)nv.size();
            nv.push_back(mesh.v[i]);
        }
    }

    for (Tri& t : mesh.f) {
        t.a = id[t.a];
        t.b = id[t.b];
        t.c = id[t.c];
    }

    mesh.v.swap(nv);
}

static bool isClosedValid(Mesh mesh) {
    compactMesh(mesh);

    if (mesh.v.empty() || mesh.f.empty()) return false;
    if (mesh.v.size() > (size_t)N) return false;

    vector<unsigned long long> edges;
    vector<unsigned long long> faces;
    edges.reserve(mesh.f.size() * 3);
    faces.reserve(mesh.f.size());

    for (const Tri& t : mesh.f) {
        if (t.a < 0 || t.b < 0 || t.c < 0) return false;
        if (t.a >= (int)mesh.v.size() || t.b >= (int)mesh.v.size() || t.c >= (int)mesh.v.size()) return false;
        if (t.a == t.b || t.a == t.c || t.b == t.c) return false;

        Vec3 cr = cross3(mesh.v[t.b] - mesh.v[t.a], mesh.v[t.c] - mesh.v[t.a]);
        if (!(norm2(cr) > areaEps2)) return false;

        faces.push_back(faceKey(t.a, t.b, t.c));
        edges.push_back(edgeKey(t.a, t.b));
        edges.push_back(edgeKey(t.b, t.c));
        edges.push_back(edgeKey(t.c, t.a));
    }

    sort(faces.begin(), faces.end());
    if (adjacent_find(faces.begin(), faces.end()) != faces.end()) return false;

    sort(edges.begin(), edges.end());
    for (size_t i = 0; i < edges.size();) {
        size_t j = i + 1;
        while (j < edges.size() && edges[j] == edges[i]) ++j;
        if (j - i != 2) return false;
        i = j;
    }

    return true;
}

struct RenderMaps {
    int r = 0;
    vector<float> depth;
    vector<float> normal;
};

static inline void projectPoint(const Vec3& p, int view, int r, double& u, double& v, double& z) {
    const double D = 2.5;
    const double f = 800.0 * (double)r / 1024.0;
    const double c = 0.5 * (double)r;

    double sx = 0.0, sy = 0.0;

    if (view == 0) {
        sx = p.y; sy = p.z; z = D - p.x;
    } else if (view == 1) {
        sx = -p.y; sy = p.z; z = D + p.x;
    } else if (view == 2) {
        sx = -p.x; sy = p.z; z = D - p.y;
    } else if (view == 3) {
        sx = p.x; sy = p.z; z = D + p.y;
    } else if (view == 4) {
        sx = p.x; sy = p.y; z = D - p.z;
    } else {
        sx = -p.x; sy = p.y; z = D + p.z;
    }

    u = f * sx / z + c;
    v = f * sy / z + c;
}

static RenderMaps renderMesh(const vector<Vec3>& v, const vector<Tri>& f, int r) {
    RenderMaps rm;
    rm.r = r;
    int pix = r * r;
    rm.depth.assign((size_t)6 * pix, 255.0f);
    rm.normal.assign((size_t)6 * pix * 3, 127.5f);

    vector<Vec3> fn(f.size());
    for (size_t i = 0; i < f.size(); ++i) {
        const Tri& t = f[i];
        fn[i] = normalized(cross3(v[t.b] - v[t.a], v[t.c] - v[t.a]));
    }

    vector<float> U(v.size()), V(v.size()), Z(v.size());

    for (int view = 0; view < 6; ++view) {
        for (size_t i = 0; i < v.size(); ++i) {
            double u, vv, z;
            projectPoint(v[i], view, r, u, vv, z);
            U[i] = (float)u;
            V[i] = (float)vv;
            Z[i] = (float)z;
        }

        float* zbuf = rm.depth.data() + (size_t)view * pix;
        float* nbuf = rm.normal.data() + (size_t)view * pix * 3;

        for (size_t ti = 0; ti < f.size(); ++ti) {
            const Tri& t = f[ti];

            int ia = t.a, ib = t.b, ic = t.c;
            float z0 = Z[ia], z1 = Z[ib], z2 = Z[ic];
            if (!(z0 > 0.0f && z1 > 0.0f && z2 > 0.0f)) continue;

            float x0 = U[ia], x1 = U[ib], x2 = U[ic];
            float y0 = V[ia], y1 = V[ib], y2 = V[ic];

            int xmin = max(0, (int)floor(min(x0, min(x1, x2))));
            int xmax = min(r - 1, (int)ceil(max(x0, max(x1, x2))));
            int ymin = max(0, (int)floor(min(y0, min(y1, y2))));
            int ymax = min(r - 1, (int)ceil(max(y0, max(y1, y2))));

            if (xmin > xmax || ymin > ymax) continue;

            float den = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2);
            if (fabs(den) < 1e-20f) continue;
            float invDen = 1.0f / den;

            Vec3 n = fn[ti];
            float nr = (float)((n.x + 1.0) * 127.5);
            float ng = (float)((n.y + 1.0) * 127.5);
            float nb = (float)((n.z + 1.0) * 127.5);

            for (int yy = ymin; yy <= ymax; ++yy) {
                float py = (float)yy + 0.5f;
                int row = yy * r;
                for (int xx = xmin; xx <= xmax; ++xx) {
                    float px = (float)xx + 0.5f;
                    float w0 = ((y1 - y2) * (px - x2) + (x2 - x1) * (py - y2)) * invDen;
                    float w1 = ((y2 - y0) * (px - x2) + (x0 - x2) * (py - y2)) * invDen;
                    float w2 = 1.0f - w0 - w1;

                    if (w0 >= -1e-6f && w1 >= -1e-6f && w2 >= -1e-6f) {
                        float depth = 1.0f / (w0 / z0 + w1 / z1 + w2 / z2);
                        int idx = row + xx;
                        if (depth < zbuf[idx]) {
                            zbuf[idx] = depth;
                            float* q = nbuf + idx * 3;
                            q[0] = nr;
                            q[1] = ng;
                            q[2] = nb;
                        }
                    }
                }
            }
        }
    }

    return rm;
}

static inline double rectSum(const vector<double>& s, int stride, int x0, int y0, int x1, int y1) {
    return s[(size_t)y1 * stride + x1]
         - s[(size_t)y0 * stride + x1]
         - s[(size_t)y1 * stride + x0]
         + s[(size_t)y0 * stride + x0];
}

static double ssimChannel(const float* x, int sx,
                          const float* y, int sy,
                          const float* dx,
                          const float* dy,
                          int r,
                          vector<double>& ix,
                          vector<double>& iy,
                          vector<double>& ix2,
                          vector<double>& iy2,
                          vector<double>& ixy) {
    int stride = r + 1;
    size_t sz = (size_t)stride * stride;

    fill(ix.begin(), ix.begin() + sz, 0.0);
    fill(iy.begin(), iy.begin() + sz, 0.0);
    fill(ix2.begin(), ix2.begin() + sz, 0.0);
    fill(iy2.begin(), iy2.begin() + sz, 0.0);
    fill(ixy.begin(), ixy.begin() + sz, 0.0);

    for (int yy = 1; yy <= r; ++yy) {
        double ax = 0, ay = 0, ax2 = 0, ay2 = 0, axy = 0;
        int row = (yy - 1) * r;

        for (int xx = 1; xx <= r; ++xx) {
            int p = row + xx - 1;
            double vx = x[(size_t)p * sx];
            double vy = y[(size_t)p * sy];

            ax += vx;
            ay += vy;
            ax2 += vx * vx;
            ay2 += vy * vy;
            axy += vx * vy;

            size_t id = (size_t)yy * stride + xx;
            size_t up = (size_t)(yy - 1) * stride + xx;

            ix[id] = ix[up] + ax;
            iy[id] = iy[up] + ay;
            ix2[id] = ix2[up] + ax2;
            iy2[id] = iy2[up] + ay2;
            ixy[id] = ixy[up] + axy;
        }
    }

    const int rad = 5;
    const int area = 121;
    const double c1 = 6.5025;
    const double c2 = 58.5225;

    long long count = 0;
    long double total = 0.0L;

    for (int yy = rad; yy < r - rad; ++yy) {
        int row = yy * r;
        for (int xx = rad; xx < r - rad; ++xx) {
            int p = row + xx;
            if (!(dx[p] < 254.0f || dy[p] < 254.0f)) continue;

            int x0 = xx - rad;
            int y0 = yy - rad;
            int x1 = xx + rad + 1;
            int y1 = yy + rad + 1;

            double sxv = rectSum(ix, stride, x0, y0, x1, y1);
            double syv = rectSum(iy, stride, x0, y0, x1, y1);
            double sx2v = rectSum(ix2, stride, x0, y0, x1, y1);
            double sy2v = rectSum(iy2, stride, x0, y0, x1, y1);
            double sxyv = rectSum(ixy, stride, x0, y0, x1, y1);

            double mux = sxv / area;
            double muy = syv / area;
            double varx = sx2v / area - mux * mux;
            double vary = sy2v / area - muy * muy;
            double cov = sxyv / area - mux * muy;

            if (varx < 0 && varx > -1e-7) varx = 0;
            if (vary < 0 && vary > -1e-7) vary = 0;

            double den = (mux * mux + muy * muy + c1) * (varx + vary + c2);
            double val = den ? ((2 * mux * muy + c1) * (2 * cov + c2) / den) : 1.0;

            total += val;
            ++count;
        }
    }

    return count ? (double)(total / count) : 1.0;
}

static double visualScore(const RenderMaps& a, const RenderMaps& b) {
    int r = a.r;
    int pix = r * r;
    int stride = r + 1;
    size_t sz = (size_t)stride * stride;

    vector<double> ix(sz), iy(sz), ix2(sz), iy2(sz), ixy(sz);

    double total = 0.0;

    for (int view = 0; view < 6; ++view) {
        const float* ad = a.depth.data() + (size_t)view * pix;
        const float* bd = b.depth.data() + (size_t)view * pix;

        double ns = 0.0;
        for (int ch = 0; ch < 3; ++ch) {
            const float* an = a.normal.data() + (size_t)view * pix * 3 + ch;
            const float* bn = b.normal.data() + (size_t)view * pix * 3 + ch;
            ns += ssimChannel(an, 3, bn, 3, ad, bd, r, ix, iy, ix2, iy2, ixy);
        }
        ns /= 3.0;

        double ds = ssimChannel(ad, 1, bd, 1, ad, bd, r, ix, iy, ix2, iy2, ixy);

        total += 0.5 * (ns + ds);
    }

    return total / 6.0;
}

static Mesh bestMesh;
static RenderMaps originalRender;
static int renderResolution = 256;

static void considerCandidate(Mesh candidate, double threshold) {
    if (elapsedSeconds() > 19.0) return;

    compactMesh(candidate);

    if (candidate.v.empty() || candidate.f.empty()) return;
    if (candidate.v.size() >= bestMesh.v.size()) return;
    if (!isClosedValid(candidate)) return;

    RenderMaps r = renderMesh(candidate.v, candidate.f, renderResolution);
    double s = visualScore(originalRender, r);

    if (s >= threshold && candidate.v.size() < bestMesh.v.size()) {
        bestMesh = move(candidate);
    }
}

static void orientByCenter(Mesh& mesh, const Vec3& center) {
    for (Tri& t : mesh.f) {
        Vec3 cr = cross3(mesh.v[t.b] - mesh.v[t.a], mesh.v[t.c] - mesh.v[t.a]);
        Vec3 c = (mesh.v[t.a] + mesh.v[t.b] + mesh.v[t.c]) / 3.0;
        if (dot3(cr, c - center) < 0) swap(t.b, t.c);
    }
}

static void jacobiAxes(Vec3 axes[3]) {
    Vec3 mean;
    for (const Vec3& p : origV) mean += p;
    mean = mean / max(1, N);

    double a[3][3] = {};
    int stride = max(1, N / 250000);
    int cnt = 0;

    for (int i = 0; i < N; i += stride) {
        Vec3 q = origV[i] - mean;
        double x[3] = {q.x, q.y, q.z};
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                a[r][c] += x[r] * x[c];
        ++cnt;
    }

    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c)
            a[r][c] /= max(1, cnt);

    double v[3][3] = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}
    };

    for (int it = 0; it < 40; ++it) {
        int p = 0, q = 1;
        double best = fabs(a[0][1]);

        if (fabs(a[0][2]) > best) {
            p = 0; q = 2; best = fabs(a[0][2]);
        }
        if (fabs(a[1][2]) > best) {
            p = 1; q = 2; best = fabs(a[1][2]);
        }

        if (best < 1e-18) break;

        double app = a[p][p];
        double aqq = a[q][q];
        double apq = a[p][q];

        double tau = (aqq - app) / (2.0 * apq);
        double t = (tau >= 0.0 ? 1.0 : -1.0) / (fabs(tau) + sqrt(1.0 + tau * tau));
        double cs = 1.0 / sqrt(1.0 + t * t);
        double sn = t * cs;

        for (int k = 0; k < 3; ++k) {
            if (k == p || k == q) continue;

            double akp = a[k][p];
            double akq = a[k][q];

            a[k][p] = a[p][k] = cs * akp - sn * akq;
            a[k][q] = a[q][k] = sn * akp + cs * akq;
        }

        a[p][p] = cs * cs * app - 2.0 * sn * cs * apq + sn * sn * aqq;
        a[q][q] = sn * sn * app + 2.0 * sn * cs * apq + cs * cs * aqq;
        a[p][q] = a[q][p] = 0.0;

        for (int k = 0; k < 3; ++k) {
            double vkp = v[k][p];
            double vkq = v[k][q];
            v[k][p] = cs * vkp - sn * vkq;
            v[k][q] = sn * vkp + cs * vkq;
        }
    }

    int order[3] = {0, 1, 2};
    sort(order, order + 3, [&](int x, int y) {
        return a[x][x] > a[y][y];
    });

    for (int j = 0; j < 3; ++j) {
        int col = order[j];
        axes[j] = normalized(Vec3(v[0][col], v[1][col], v[2][col]));
    }

    if (dot3(cross3(axes[0], axes[1]), axes[2]) < 0) axes[2] = axes[2] * -1.0;
}

static bool ellipsoidFit(const Vec3 axes[3], Vec3& center, double radius[3], double& rms, double& maxErr) {
    double lo[3] = {1e100, 1e100, 1e100};
    double hi[3] = {-1e100, -1e100, -1e100};

    for (const Vec3& p : origV) {
        for (int k = 0; k < 3; ++k) {
            double t = dot3(p, axes[k]);
            lo[k] = min(lo[k], t);
            hi[k] = max(hi[k], t);
        }
    }

    center = Vec3();

    for (int k = 0; k < 3; ++k) {
        double mid = 0.5 * (lo[k] + hi[k]);
        radius[k] = 0.5 * (hi[k] - lo[k]);
        if (!(radius[k] > 1e-10 * bboxDiag)) return false;
        center += axes[k] * mid;
    }

    int stride = max(1, N / 300000);
    int cnt = 0;
    double sumSq = 0.0;
    maxErr = 0.0;

    for (int i = 0; i < N; i += stride) {
        Vec3 q = origV[i] - center;
        double r2 = 0.0;

        for (int k = 0; k < 3; ++k) {
            double u = dot3(q, axes[k]) / radius[k];
            r2 += u * u;
        }

        double e = fabs(sqrt(max(0.0, r2)) - 1.0);
        sumSq += e * e;
        maxErr = max(maxErr, e);
        ++cnt;
    }

    if (cnt == 0) return false;

    rms = sqrt(sumSq / cnt);
    return true;
}

static Mesh makeEllipsoidMesh(const Vec3 axes[3],
                              const Vec3& center,
                              const double radius[3],
                              int lat,
                              int lon) {
    Mesh mesh;
    mesh.v.reserve(2 + (lat - 1) * lon);
    mesh.f.reserve(2 * lat * lon);

    const double pi = acos(-1.0);

    auto point = [&](double x, double y, double z) {
        return center + axes[0] * (radius[0] * x)
                      + axes[1] * (radius[1] * y)
                      + axes[2] * (radius[2] * z);
    };

    mesh.v.push_back(point(0, 0, 1));
    mesh.v.push_back(point(0, 0, -1));

    auto id = [&](int r, int j) {
        return 2 + (r - 1) * lon + (j % lon + lon) % lon;
    };

    for (int r = 1; r < lat; ++r) {
        double theta = pi * r / lat;
        double st = sin(theta);
        double ct = cos(theta);

        for (int j = 0; j < lon; ++j) {
            double phi = 2.0 * pi * j / lon;
            mesh.v.push_back(point(st * cos(phi), st * sin(phi), ct));
        }
    }

    for (int j = 0; j < lon; ++j) {
        mesh.f.push_back({0, id(1, j), id(1, j + 1)});
    }

    for (int r = 1; r < lat - 1; ++r) {
        for (int j = 0; j < lon; ++j) {
            int a = id(r, j);
            int b = id(r + 1, j);
            int c = id(r + 1, j + 1);
            int d = id(r, j + 1);

            mesh.f.push_back({a, b, c});
            mesh.f.push_back({a, c, d});
        }
    }

    for (int j = 0; j < lon; ++j) {
        mesh.f.push_back({1, id(lat - 1, j + 1), id(lat - 1, j)});
    }

    orientByCenter(mesh, center);
    return mesh;
}

static void tryEllipsoidCandidates() {
    if (N < 900 || elapsedSeconds() > 5.0) return;

    Vec3 identity[3] = {Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1)};
    Vec3 pca[3];
    jacobiAxes(pca);

    for (int type = 0; type < 2; ++type) {
        Vec3* axes = type == 0 ? identity : pca;

        Vec3 center;
        double radius[3];
        double rms, maxErr;

        if (!ellipsoidFit(axes, center, radius, rms, maxErr)) continue;

        double rmax = max(radius[0], max(radius[1], radius[2]));
        double rmin = min(radius[0], min(radius[1], radius[2]));

        if (rmin < 0.14 * rmax) continue;

        double rmsLimit = N < 5000 ? 0.0065 : 0.0048;
        double maxLimit = N < 5000 ? 0.028 : 0.020;

        if (rms > rmsLimit || maxErr > maxLimit) continue;

        vector<pair<int, int>> trials = {
            {12, 24}, {16, 32}, {18, 36}, {20, 40}, {24, 48}, {28, 56}
        };

        for (auto pr : trials) {
            if (elapsedSeconds() > 10.0) break;

            Mesh mesh = makeEllipsoidMesh(axes, center, radius, pr.first, pr.second);
            if (mesh.v.size() < bestMesh.v.size()) {
                considerCandidate(move(mesh), 0.920);
            }
        }
    }
}

static vector<Vec3> originalVertexNormals() {
    vector<Vec3> n(N);
    for (const Tri& t : origF) {
        Vec3 cr = cross3(origV[t.b] - origV[t.a], origV[t.c] - origV[t.a]);
        n[t.a] += cr;
        n[t.b] += cr;
        n[t.c] += cr;
    }
    return n;
}

static bool faceUsesOnly(const Tri& f, int a, int b, int c, int d) {
    auto ok = [&](int x) {
        return x == a || x == b || x == c || x == d;
    };
    return ok(f.a) && ok(f.b) && ok(f.c);
}

static void addFaceOriented(Mesh& mesh, int a, int b, int c, const Vec3& ref) {
    Vec3 cr = cross3(mesh.v[b] - mesh.v[a], mesh.v[c] - mesh.v[a]);
    if (dot3(cr, ref) < 0) swap(b, c);
    mesh.f.push_back({a, b, c});
}

static void tryOrderedTorusGrid() {
    if (N < 300 || M != 2 * N || elapsedSeconds() > 10.0) return;

    vector<int> candidates;
    for (int s = 6; s * s <= N; ++s) {
        if (N % s == 0) {
            candidates.push_back(s);
            candidates.push_back(N / s);
        }
    }

    sort(candidates.begin(), candidates.end());
    candidates.erase(unique(candidates.begin(), candidates.end()), candidates.end());

    int S = 0, U = 0;

    for (int s : candidates) {
        int u = N / s;
        if (s < 6 || u < 6 || s > 5000 || u > 5000) continue;

        int step = max(1, N / 1200);
        int total = 0, ok = 0;

        for (int cell = 0; cell < N; cell += step) {
            int i = cell / s;
            int j = cell - i * s;

            int a = i * s + j;
            int b = i * s + (j + 1) % s;
            int c = ((i + 1) % u) * s + j;
            int d = ((i + 1) % u) * s + (j + 1) % s;

            int f0 = 2 * cell;
            int f1 = f0 + 1;

            if (f1 < M && faceUsesOnly(origF[f0], a, b, c, d) && faceUsesOnly(origF[f1], a, b, c, d)) ++ok;
            ++total;
        }

        if (total > 100 && ok * 1000 >= total * 990) {
            S = s;
            U = u;
            break;
        }
    }

    if (!S) return;

    vector<Vec3> vn = originalVertexNormals();
    vector<int> targets = {768, 1024, 1536, 2048, 3072, 4096, 6144, 8192, 12288, 16384};

    for (int target : targets) {
        if (elapsedSeconds() > 15.5) break;

        double aspect = sqrt((double)U / (double)S);
        int U2 = max(6, min(U, (int)(sqrt((double)target) * aspect + 0.5)));
        int S2 = max(6, min(S, target / max(1, U2)));

        if (U2 * S2 >= (int)bestMesh.v.size() || U2 * S2 >= N) continue;

        Mesh mesh;
        mesh.v.reserve(U2 * S2);
        mesh.f.reserve(2 * U2 * S2);

        vector<int> src;
        src.reserve(U2 * S2);

        for (int i = 0; i < U2; ++i) {
            int oi = (int)((long long)i * U / U2);
            for (int j = 0; j < S2; ++j) {
                int oj = (int)((long long)j * S / S2);
                int old = oi * S + oj;
                src.push_back(old);
                mesh.v.push_back(origV[old]);
            }
        }

        auto id = [&](int i, int j) {
            i = (i % U2 + U2) % U2;
            j = (j % S2 + S2) % S2;
            return i * S2 + j;
        };

        for (int i = 0; i < U2; ++i) {
            for (int j = 0; j < S2; ++j) {
                int a = id(i, j);
                int b = id(i + 1, j);
                int c = id(i + 1, j + 1);
                int d = id(i, j + 1);

                Vec3 ref1 = vn[src[a]] + vn[src[b]] + vn[src[d]];
                Vec3 ref2 = vn[src[b]] + vn[src[c]] + vn[src[d]];

                addFaceOriented(mesh, a, b, d, ref1);
                addFaceOriented(mesh, b, c, d, ref2);
            }
        }

        considerCandidate(move(mesh), 0.920);
    }
}

static bool sameTriangleSet(const Tri& f, int a, int b, int c) {
    int x[3] = {f.a, f.b, f.c};
    int y[3] = {a, b, c};
    sort(x, x + 3);
    sort(y, y + 3);
    return x[0] == y[0] && x[1] == y[1] && x[2] == y[2];
}

static void tryOrderedSphereGrid() {
    if (N < 300 || M != 2 * (N - 2) || elapsedSeconds() > 12.0) return;

    int S = 0, R = 0;

    for (int s = 8; s <= min(4096, N - 2); ++s) {
        if ((N - 2) % s != 0) continue;

        int r = (N - 2) / s;
        if (r < 3) continue;

        int step = max(1, (r - 1) * s / 1000);
        int total = 0, ok = 0;

        for (int q = 0; q < (r - 1) * s; q += step) {
            int rr = q / s;
            int j = q - rr * s;

            if (rr == 0) {
                int fid = j;
                int a = 1 + j;
                int b = 1 + (j + 1) % s;
                if (fid < M && (origF[fid].a == 0 || origF[fid].b == 0 || origF[fid].c == 0)
                    && sameTriangleSet(origF[fid], 0, a, b)) {
                    ++ok;
                }
            } else if (rr < r - 1) {
                int fid = s + 2 * ((rr - 1) * s + j);

                int a = 1 + (rr - 1) * s + j;
                int b = 1 + (rr - 1) * s + (j + 1) % s;
                int c = 1 + rr * s + j;
                int d = 1 + rr * s + (j + 1) % s;

                if (fid + 1 < M && faceUsesOnly(origF[fid], a, b, c, d)
                    && faceUsesOnly(origF[fid + 1], a, b, c, d)) {
                    ++ok;
                }
            }

            ++total;
        }

        if (total > 100 && ok * 100 >= total * 95) {
            S = s;
            R = r;
            break;
        }
    }

    if (!S) return;

    vector<Vec3> vn = originalVertexNormals();

    vector<pair<int, int>> trials = {
        {max(4, R / 10), max(12, S / 10)},
        {max(5, R / 8), max(16, S / 8)},
        {max(6, R / 6), max(20, S / 6)},
        {max(8, R / 4), max(24, S / 4)}
    };

    for (auto pr : trials) {
        if (elapsedSeconds() > 16.0) break;

        int R2 = pr.first;
        int S2 = pr.second;

        if (R2 < 3 || S2 < 8) continue;
        if (2 + (R2 - 1) * S2 >= (int)bestMesh.v.size()) continue;
        if (2 + (R2 - 1) * S2 >= N) continue;

        Mesh mesh;
        vector<int> src;

        mesh.v.reserve(2 + (R2 - 1) * S2);
        mesh.f.reserve(2 * R2 * S2);

        mesh.v.push_back(origV[0]);
        src.push_back(0);

        for (int r = 1; r < R2; ++r) {
            int orr = 1 + (int)((long long)(r - 1) * (R - 1) / max(1, R2 - 1));

            for (int j = 0; j < S2; ++j) {
                int oj = (int)((long long)j * S / S2);
                int old = 1 + (orr - 1) * S + oj;
                src.push_back(old);
                mesh.v.push_back(origV[old]);
            }
        }

        int bottom = (int)mesh.v.size();
        mesh.v.push_back(origV[N - 1]);
        src.push_back(N - 1);

        auto id = [&](int r, int j) {
            return 1 + (r - 1) * S2 + (j % S2 + S2) % S2;
        };

        for (int j = 0; j < S2; ++j) {
            addFaceOriented(mesh, 0, id(1, j + 1), id(1, j), vn[0]);
        }

        for (int r = 1; r < R2 - 1; ++r) {
            for (int j = 0; j < S2; ++j) {
                int a = id(r, j);
                int b = id(r, j + 1);
                int c = id(r + 1, j);
                int d = id(r + 1, j + 1);

                addFaceOriented(mesh, a, b, c, vn[src[a]] + vn[src[b]] + vn[src[c]]);
                addFaceOriented(mesh, b, d, c, vn[src[b]] + vn[src[d]] + vn[src[c]]);
            }
        }

        for (int j = 0; j < S2; ++j) {
            addFaceOriented(mesh, bottom, id(R2 - 1, j), id(R2 - 1, j + 1), vn[N - 1]);
        }

        considerCandidate(move(mesh), 0.920);
    }
}

struct CollapseVertex {
    Vec3 p;
    Quadric q;
    double radius = 0.0;
    int version = 0;
    bool alive = true;
    vector<int> inc;
};

struct CollapseFace {
    int a, b, c;
    Vec3 n;
    bool alive = true;
};

struct QueueNode {
    double cost;
    int a, b;
    int va, vb;
    bool operator<(const QueueNode& o) const {
        return cost > o.cost;
    }
};

static vector<CollapseVertex> cv;
static vector<CollapseFace> cf;
static int aliveVertices, aliveFaces;
static vector<int> marker;
static int markerToken = 1;

static inline bool faceHas(const CollapseFace& f, int v) {
    return f.a == v || f.b == v || f.c == v;
}

static Vec3 faceNormalCurrent(int fid) {
    const CollapseFace& f = cf[fid];
    return normalized(cross3(cv[f.b].p - cv[f.a].p, cv[f.c].p - cv[f.a].p));
}

static void cleanupIncident(int v) {
    if (!cv[v].alive) {
        vector<int>().swap(cv[v].inc);
        return;
    }

    if (cv[v].inc.size() < 128) return;

    vector<int> keep;
    keep.reserve(cv[v].inc.size());

    for (int fid : cv[v].inc) {
        if (fid >= 0 && fid < (int)cf.size() && cf[fid].alive && faceHas(cf[fid], v)) {
            keep.push_back(fid);
        }
    }

    cv[v].inc.swap(keep);
}

static bool solveOptimal(const Quadric& q, Vec3& out) {
    double a00 = q.q[0], a01 = q.q[1], a02 = q.q[2];
    double a11 = q.q[4], a12 = q.q[5], a22 = q.q[7];
    double b0 = -q.q[3], b1 = -q.q[6], b2 = -q.q[8];

    double det = a00 * (a11 * a22 - a12 * a12)
               - a01 * (a01 * a22 - a12 * a02)
               + a02 * (a01 * a12 - a11 * a02);

    if (fabs(det) < 1e-14) return false;

    double dx = b0 * (a11 * a22 - a12 * a12)
              - a01 * (b1 * a22 - a12 * b2)
              + a02 * (b1 * a12 - a11 * b2);

    double dy = a00 * (b1 * a22 - a12 * b2)
              - b0 * (a01 * a22 - a12 * a02)
              + a02 * (a01 * b2 - b1 * a02);

    double dz = a00 * (a11 * b2 - b1 * a12)
              - a01 * (a01 * b2 - b1 * a02)
              + b0 * (a01 * a12 - a11 * a02);

    out = Vec3(dx / det, dy / det, dz / det);

    return isfinite(out.x) && isfinite(out.y) && isfinite(out.z);
}

static bool collectOpposites(int a, int b, int opp[2]) {
    cleanupIncident(a);
    int cnt = 0;

    for (int fid : cv[a].inc) {
        if (!cf[fid].alive) continue;

        const CollapseFace& f = cf[fid];
        if (!faceHas(f, a) || !faceHas(f, b)) continue;

        if (cnt >= 2) return false;

        int t = -1;
        if (f.a != a && f.a != b) t = f.a;
        if (f.b != a && f.b != b) t = f.b;
        if (f.c != a && f.c != b) t = f.c;

        if (t < 0) return false;
        opp[cnt++] = t;
    }

    return cnt == 2 && opp[0] != opp[1];
}

static bool linkOK(int a, int b) {
    int opp[2] = {-1, -1};
    if (!collectOpposites(a, b, opp)) return false;

    if (++markerToken > 1000000000) {
        fill(marker.begin(), marker.end(), 0);
        markerToken = 1;
    }

    int tokenA = markerToken++;
    int tokenB = markerToken++;

    for (int fid : cv[a].inc) {
        if (!cf[fid].alive || !faceHas(cf[fid], a)) continue;

        const CollapseFace& f = cf[fid];
        int vv[3] = {f.a, f.b, f.c};

        for (int x : vv) {
            if (x != a && x != b) marker[x] = tokenA;
        }
    }

    int common = 0;
    int got0 = 0, got1 = 0;

    cleanupIncident(b);

    for (int fid : cv[b].inc) {
        if (!cf[fid].alive || !faceHas(cf[fid], b)) continue;

        const CollapseFace& f = cf[fid];
        int vv[3] = {f.a, f.b, f.c};

        for (int x : vv) {
            if (x == a || x == b) continue;

            if (marker[x] == tokenA) {
                marker[x] = tokenB;
                ++common;
                if (x == opp[0]) got0 = 1;
                if (x == opp[1]) got1 = 1;
                if (common > 2) return false;
            }
        }
    }

    return common == 2 && got0 && got1;
}

static bool flipOK(int a, int b, const Vec3& pos) {
    static vector<int> touched;
    touched.clear();

    for (int fid : cv[a].inc) {
        if (cf[fid].alive) touched.push_back(fid);
    }
    for (int fid : cv[b].inc) {
        if (cf[fid].alive) touched.push_back(fid);
    }

    sort(touched.begin(), touched.end());
    touched.erase(unique(touched.begin(), touched.end()), touched.end());

    double minCos = cos(28.0 * acos(-1.0) / 180.0);

    for (int fid : touched) {
        CollapseFace& f = cf[fid];

        bool hasA = faceHas(f, a);
        bool hasB = faceHas(f, b);

        if (hasA && hasB) continue;

        Vec3 pa = (f.a == a || f.a == b) ? pos : cv[f.a].p;
        Vec3 pb = (f.b == a || f.b == b) ? pos : cv[f.b].p;
        Vec3 pc = (f.c == a || f.c == b) ? pos : cv[f.c].p;

        Vec3 cr = cross3(pb - pa, pc - pa);
        double l = norm3(cr);

        if (!(l > 1e-18 * bboxDiag * bboxDiag)) return false;

        Vec3 n = cr / l;
        if (dot3(n, f.n) < minCos) return false;
    }

    return true;
}

static bool bestCollapsePosition(int a, int b, Vec3& pos, double& rad, double& cost) {
    if (a == b || !cv[a].alive || !cv[b].alive) return false;
    if (!linkOK(a, b)) return false;

    Quadric q = cv[a].q;
    q.add(cv[b].q);

    Vec3 pa = cv[a].p;
    Vec3 pb = cv[b].p;

    Vec3 candidates[7];
    int cnt = 0;

    Vec3 opt;
    if (solveOptimal(q, opt)) candidates[cnt++] = opt;

    candidates[cnt++] = (pa + pb) * 0.5;
    candidates[cnt++] = pa;
    candidates[cnt++] = pb;
    candidates[cnt++] = pa * 0.75 + pb * 0.25;
    candidates[cnt++] = pa * 0.25 + pb * 0.75;

    double limit = 0.0492 * bboxDiag;

    bool ok = false;
    cost = 1e100;

    for (int i = 0; i < cnt; ++i) {
        Vec3 p = candidates[i];

        if (!isfinite(p.x) || !isfinite(p.y) || !isfinite(p.z)) continue;

        double r = max(cv[a].radius + norm3(p - pa), cv[b].radius + norm3(p - pb));
        if (r > limit) continue;

        if (!flipOK(a, b, p)) continue;

        double c = q.eval(p) + 1e-10 * norm2(pa - pb);
        if (c < cost) {
            cost = c;
            pos = p;
            rad = r;
            ok = true;
        }
    }

    return ok;
}

static double cheapCost(int a, int b) {
    Quadric q = cv[a].q;
    q.add(cv[b].q);

    Vec3 mid = (cv[a].p + cv[b].p) * 0.5;
    Vec3 opt;

    double c = q.eval(mid);
    if (solveOptimal(q, opt)) c = min(c, q.eval(opt));

    return c + 1e-10 * norm2(cv[a].p - cv[b].p);
}

static void pushEdge(priority_queue<QueueNode>& pq, int a, int b) {
    if (a == b || !cv[a].alive || !cv[b].alive) return;
    if (a > b) swap(a, b);

    pq.push({cheapCost(a, b), a, b, cv[a].version, cv[b].version});
}

static void performCollapse(int keep, int remove, const Vec3& pos, double radius, priority_queue<QueueNode>& pq) {
    cleanupIncident(keep);
    cleanupIncident(remove);

    cv[keep].q.add(cv[remove].q);
    cv[keep].p = pos;
    cv[keep].radius = radius;

    vector<int> incidentRemove = cv[remove].inc;

    for (int fid : incidentRemove) {
        if (!cf[fid].alive) continue;

        CollapseFace& f = cf[fid];
        if (!faceHas(f, remove)) continue;

        bool hasKeep = faceHas(f, keep);

        if (hasKeep) {
            f.alive = false;
            --aliveFaces;
        } else {
            if (f.a == remove) f.a = keep;
            if (f.b == remove) f.b = keep;
            if (f.c == remove) f.c = keep;

            f.n = faceNormalCurrent(fid);
            cv[keep].inc.push_back(fid);
        }
    }

    cv[remove].alive = false;
    vector<int>().swap(cv[remove].inc);

    --aliveVertices;
    ++cv[keep].version;
    ++cv[remove].version;

    for (int fid : cv[keep].inc) {
        if (cf[fid].alive && faceHas(cf[fid], keep)) {
            cf[fid].n = faceNormalCurrent(fid);
        }
    }

    cleanupIncident(keep);

    if (++markerToken > 1000000000) {
        fill(marker.begin(), marker.end(), 0);
        markerToken = 1;
    }

    int tok = markerToken++;
    vector<int> neigh;

    for (int fid : cv[keep].inc) {
        if (!cf[fid].alive) continue;

        int vv[3] = {cf[fid].a, cf[fid].b, cf[fid].c};

        for (int x : vv) {
            if (x != keep && cv[x].alive && marker[x] != tok) {
                marker[x] = tok;
                neigh.push_back(x);
            }
        }
    }

    for (int x : neigh) pushEdge(pq, keep, x);
}

static Mesh snapshotCurrentCollapseMesh() {
    Mesh mesh;
    vector<int> id(cv.size(), -1);

    for (int i = 0; i < (int)cv.size(); ++i) {
        if (cv[i].alive) {
            id[i] = (int)mesh.v.size();
            mesh.v.push_back(cv[i].p);
        }
    }

    for (const CollapseFace& f : cf) {
        if (!f.alive) continue;
        if (f.a == f.b || f.a == f.c || f.b == f.c) continue;
        if (id[f.a] < 0 || id[f.b] < 0 || id[f.c] < 0) continue;

        mesh.f.push_back({id[f.a], id[f.b], id[f.c]});
    }

    compactMesh(mesh);
    return mesh;
}

static void tryConservativeQEM() {
    if (N < 50 || elapsedSeconds() > 13.0) return;

    cv.assign(N, CollapseVertex());
    cf.assign(M, CollapseFace());

    vector<int> deg(N, 0);

    for (int i = 0; i < N; ++i) {
        cv[i].p = origV[i];
    }

    for (int i = 0; i < M; ++i) {
        cf[i].a = origF[i].a;
        cf[i].b = origF[i].b;
        cf[i].c = origF[i].c;
        cf[i].alive = true;

        ++deg[cf[i].a];
        ++deg[cf[i].b];
        ++deg[cf[i].c];
    }

    for (int i = 0; i < N; ++i) cv[i].inc.reserve(deg[i] + 8);

    for (int i = 0; i < M; ++i) {
        cv[cf[i].a].inc.push_back(i);
        cv[cf[i].b].inc.push_back(i);
        cv[cf[i].c].inc.push_back(i);
        cf[i].n = faceNormalCurrent(i);
    }

    for (int i = 0; i < M; ++i) {
        CollapseFace& f = cf[i];

        Vec3 a = cv[f.a].p;
        Vec3 b = cv[f.b].p;
        Vec3 c = cv[f.c].p;

        Vec3 cr = cross3(b - a, c - a);
        double l = norm3(cr);
        if (!(l > 0.0)) continue;

        Vec3 n = cr / l;
        double d = -dot3(n, a);
        double w = max(1e-12, 0.5 * l);

        cv[f.a].q.addPlane(n.x, n.y, n.z, d, w);
        cv[f.b].q.addPlane(n.x, n.y, n.z, d, w);
        cv[f.c].q.addPlane(n.x, n.y, n.z, d, w);
    }

    aliveVertices = N;
    aliveFaces = M;

    marker.assign(N, 0);
    markerToken = 1;

    priority_queue<QueueNode> pq;

    vector<unsigned long long> edges;
    edges.reserve((size_t)M * 3);

    for (const Tri& t : origF) {
        edges.push_back(edgeKey(t.a, t.b));
        edges.push_back(edgeKey(t.b, t.c));
        edges.push_back(edgeKey(t.c, t.a));
    }

    sort(edges.begin(), edges.end());
    edges.erase(unique(edges.begin(), edges.end()), edges.end());

    for (unsigned long long k : edges) {
        pushEdge(pq, (int)(k >> 32), (int)(k & 0xffffffffu));
    }

    vector<int> targets;
    for (double ratio : {0.55, 0.42, 0.32, 0.25, 0.20, 0.16, 0.13, 0.11, 0.095, 0.08, 0.065, 0.05}) {
        int t = max(4, (int)ceil((double)N * ratio));
        if (t < N) targets.push_back(t);
    }

    sort(targets.begin(), targets.end(), greater<int>());

    int targetIndex = 0;
    vector<Mesh> snapshots;

    while (!pq.empty() && targetIndex < (int)targets.size() && elapsedSeconds() < 16.2) {
        QueueNode node = pq.top();
        pq.pop();

        int a = node.a;
        int b = node.b;

        if (a == b || !cv[a].alive || !cv[b].alive) continue;
        if (node.va != cv[a].version || node.vb != cv[b].version) continue;

        Vec3 pos;
        double rad = 0.0, cost = 0.0;

        if (!bestCollapsePosition(a, b, pos, rad, cost)) continue;

        int keep = a;
        int rem = b;

        if (cv[b].inc.size() > cv[a].inc.size()) {
            keep = b;
            rem = a;
        }

        performCollapse(keep, rem, pos, rad, pq);

        while (targetIndex < (int)targets.size() && aliveVertices <= targets[targetIndex]) {
            snapshots.push_back(snapshotCurrentCollapseMesh());
            ++targetIndex;
        }
    }

    if (snapshots.empty()) snapshots.push_back(snapshotCurrentCollapseMesh());

    sort(snapshots.begin(), snapshots.end(), [](const Mesh& a, const Mesh& b) {
        return a.v.size() < b.v.size();
    });

    for (Mesh& mesh : snapshots) {
        if (elapsedSeconds() > 19.0) break;
        if (mesh.v.size() < bestMesh.v.size()) {
            considerCandidate(move(mesh), 0.918);
        }
    }
}

static double pointTriangleDistance2(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c) {
    Vec3 ab = b - a;
    Vec3 ac = c - a;
    Vec3 ap = p - a;

    double d1 = dot3(ab, ap);
    double d2 = dot3(ac, ap);

    if (d1 <= 0.0 && d2 <= 0.0) return norm2(ap);

    Vec3 bp = p - b;
    double d3 = dot3(ab, bp);
    double d4 = dot3(ac, bp);

    if (d3 >= 0.0 && d4 <= d3) return norm2(bp);

    double vc = d1 * d4 - d3 * d2;

    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        double v = d1 / (d1 - d3);
        Vec3 q = a + ab * v;
        return norm2(p - q);
    }

    Vec3 cp = p - c;
    double d5 = dot3(ab, cp);
    double d6 = dot3(ac, cp);

    if (d6 >= 0.0 && d5 <= d6) return norm2(cp);

    double vb = d5 * d2 - d1 * d6;

    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        double w = d2 / (d2 - d6);
        Vec3 q = a + ac * w;
        return norm2(p - q);
    }

    double va = d3 * d6 - d5 * d4;

    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        Vec3 q = b + (c - b) * w;
        return norm2(p - q);
    }

    Vec3 n = cross3(ab, ac);
    double nn = norm2(n);
    if (!(nn > 0.0)) return 1e100;

    double dist = dot3(p - a, n);
    return dist * dist / nn;
}

static void tryLocalVertexRemoval() {
    if (bestMesh.v.size() >= origV.size() || bestMesh.v.size() < 20 || elapsedSeconds() > 17.0) return;

    Mesh cur = bestMesh;

    int n = (int)cur.v.size();
    int f0 = (int)cur.f.size();

    vector<char> aliveV(n, 1);
    vector<char> aliveF(f0, 1);
    vector<vector<int>> inc(n);

    for (int i = 0; i < f0; ++i) {
        inc[cur.f[i].a].push_back(i);
        inc[cur.f[i].b].push_back(i);
        inc[cur.f[i].c].push_back(i);
    }

    auto has = [&](int fid, int v) {
        const Tri& t = cur.f[fid];
        return t.a == v || t.b == v || t.c == v;
    };

    struct Patch {
        int v = -1;
        double score = 1e100;
        vector<int> ring;
        vector<int> oldFaces;
        vector<Tri> newFaces;
    };

    auto inspect = [&](int v, Patch& patch) -> bool {
        patch = Patch();
        patch.v = v;

        if (v < 0 || v >= n || !aliveV[v]) return false;

        vector<int> old;
        for (int fid : inc[v]) {
            if (fid >= 0 && fid < (int)cur.f.size() && aliveF[fid] && has(fid, v)) {
                old.push_back(fid);
            }
        }

        sort(old.begin(), old.end());
        old.erase(unique(old.begin(), old.end()), old.end());

        int k = (int)old.size();
        if (k < 3 || k > 8) return false;

        vector<int> nb;
        vector<pair<int, int>> boundary;
        Vec3 normalSum;
        double maxRadius = 0.0;

        for (int fid : old) {
            const Tri& t = cur.f[fid];
            int all[3] = {t.a, t.b, t.c};
            int other[2];
            int cnt = 0;

            for (int x : all) {
                if (x != v) other[cnt++] = x;
            }

            if (cnt != 2 || other[0] == other[1] || !aliveV[other[0]] || !aliveV[other[1]]) return false;

            for (int z : other) {
                if (find(nb.begin(), nb.end(), z) == nb.end()) nb.push_back(z);
            }

            boundary.push_back({other[0], other[1]});

            Vec3 cr = cross3(cur.v[t.b] - cur.v[t.a], cur.v[t.c] - cur.v[t.a]);
            double l = norm3(cr);
            if (!(l > 0.0)) return false;
            normalSum += cr / l;

            maxRadius = max(maxRadius, norm3(cur.v[other[0]] - cur.v[v]));
            maxRadius = max(maxRadius, norm3(cur.v[other[1]] - cur.v[v]));
        }

        if ((int)nb.size() != k) return false;
        if (maxRadius > 0.090 * bboxDiag) return false;

        double nl = norm3(normalSum);
        if (!(nl > 1e-12)) return false;

        Vec3 normal = normalSum / nl;

        double oldCos = cos(17.0 * acos(-1.0) / 180.0);
        for (int fid : old) {
            const Tri& t = cur.f[fid];
            Vec3 cr = cross3(cur.v[t.b] - cur.v[t.a], cur.v[t.c] - cur.v[t.a]);
            double l = norm3(cr);
            if (l <= 0.0 || dot3(cr / l, normal) < oldCos) return false;
        }

        vector<vector<int>> adj(k);
        auto indexOf = [&](int x) {
            for (int i = 0; i < k; ++i) if (nb[i] == x) return i;
            return -1;
        };

        set<pair<int, int>> boundarySet;

        for (auto e : boundary) {
            int a = indexOf(e.first);
            int b = indexOf(e.second);

            if (a < 0 || b < 0 || a == b) return false;
            if (a > b) swap(a, b);

            if (boundarySet.count({a, b})) return false;
            boundarySet.insert({a, b});

            adj[a].push_back(b);
            adj[b].push_back(a);
        }

        for (int i = 0; i < k; ++i) {
            if (adj[i].size() != 2) return false;
        }

        vector<int> order;
        order.reserve(k);

        int prev = -1;
        int curNode = 0;

        for (int step = 0; step < k; ++step) {
            order.push_back(curNode);
            int next = adj[curNode][0] == prev ? adj[curNode][1] : adj[curNode][0];
            prev = curNode;
            curNode = next;
        }

        if (curNode != 0) return false;

        for (int id : order) patch.ring.push_back(nb[id]);

        unordered_set<int> skip(old.begin(), old.end());

        auto edgeCountOutside = [&](int a, int b) {
            int c = 0;
            for (int fid : inc[a]) {
                if (fid < 0 || fid >= (int)cur.f.size() || !aliveF[fid] || skip.count(fid)) continue;
                if (has(fid, a) && has(fid, b)) ++c;
            }
            return c;
        };

        auto faceExistsOutside = [&](int a, int b, int c) {
            int x[3] = {a, b, c};
            sort(x, x + 3);

            int scan = a;
            if (inc[b].size() < inc[scan].size()) scan = b;
            if (inc[c].size() < inc[scan].size()) scan = c;

            for (int fid : inc[scan]) {
                if (fid < 0 || fid >= (int)cur.f.size() || !aliveF[fid] || skip.count(fid)) continue;

                Tri t = cur.f[fid];
                int y[3] = {t.a, t.b, t.c};
                sort(y, y + 3);

                if (x[0] == y[0] && x[1] == y[1] && x[2] == y[2]) return true;
            }

            return false;
        };

        patch.newFaces.clear();

        double bestDist2 = 1e100;
        double newCos = cos(26.0 * acos(-1.0) / 180.0);

        for (int i = 1; i + 1 < k; ++i) {
            Tri t{patch.ring[0], patch.ring[i], patch.ring[i + 1]};

            Vec3 cr = cross3(cur.v[t.b] - cur.v[t.a], cur.v[t.c] - cur.v[t.a]);

            if (dot3(cr, normal) < 0) {
                swap(t.b, t.c);
                cr = cr * -1.0;
            }

            double l = norm3(cr);
            if (!(l > 1e-12 * bboxDiag * bboxDiag)) return false;
            if (dot3(cr / l, normal) < newCos) return false;

            if (faceExistsOutside(t.a, t.b, t.c)) return false;

            int edges[3][2] = {
                {t.a, t.b},
                {t.b, t.c},
                {t.c, t.a}
            };

            for (auto& e : edges) {
                int ia = indexOf(e[0]);
                int ib = indexOf(e[1]);

                if (ia < 0 || ib < 0) return false;
                if (ia > ib) swap(ia, ib);

                bool boundaryEdge = boundarySet.count({ia, ib});
                int cnt = edgeCountOutside(e[0], e[1]);

                if (boundaryEdge) {
                    if (cnt != 1) return false;
                } else {
                    if (cnt != 0) return false;
                }
            }

            bestDist2 = min(bestDist2, pointTriangleDistance2(cur.v[v], cur.v[t.a], cur.v[t.b], cur.v[t.c]));
            patch.newFaces.push_back(t);
        }

        double d = sqrt(max(0.0, bestDist2));
        if (d > 0.043 * bboxDiag) return false;

        patch.oldFaces = old;
        patch.score = d / bboxDiag + 0.001 * k;
        return true;
    };

    vector<Patch> candidates;
    candidates.reserve(min(n, 50000));

    for (int v = 0; v < n; ++v) {
        if ((v & 1023) == 0 && elapsedSeconds() > 18.0) break;

        Patch p;
        if (inspect(v, p)) candidates.push_back(move(p));
    }

    if (candidates.empty()) return;

    sort(candidates.begin(), candidates.end(), [](const Patch& a, const Patch& b) {
        if (a.score != b.score) return a.score < b.score;
        return a.v < b.v;
    });

    int removed = 0;
    int cap = max(10, min(3000, n / 25));

    for (const Patch& p : candidates) {
        if (removed >= cap || elapsedSeconds() > 18.6) break;

        Patch fresh;
        if (!inspect(p.v, fresh)) continue;
        if (fresh.ring != p.ring) continue;

        for (int fid : fresh.oldFaces) {
            if (aliveF[fid]) aliveF[fid] = 0;
        }

        aliveV[p.v] = 0;

        for (Tri t : fresh.newFaces) {
            int id = (int)cur.f.size();
            cur.f.push_back(t);
            aliveF.push_back(1);

            inc[t.a].push_back(id);
            inc[t.b].push_back(id);
            inc[t.c].push_back(id);
        }

        ++removed;
    }

    if (!removed) return;

    Mesh out;
    vector<int> id(cur.v.size(), -1);

    for (int i = 0; i < (int)cur.v.size(); ++i) {
        if (aliveV[i]) {
            id[i] = (int)out.v.size();
            out.v.push_back(cur.v[i]);
        }
    }

    for (int i = 0; i < (int)cur.f.size(); ++i) {
        if (!aliveF[i]) continue;

        Tri t = cur.f[i];
        if (id[t.a] >= 0 && id[t.b] >= 0 && id[t.c] >= 0) {
            out.f.push_back({id[t.a], id[t.b], id[t.c]});
        }
    }

    considerCandidate(move(out), 0.917);
}

static void outputMesh(Mesh mesh) {
    compactMesh(mesh);

    if (!isClosedValid(mesh)) {
        mesh.v = origV;
        mesh.f = origF;
    }

    static char outBuf[1 << 20];
    setvbuf(stdout, outBuf, _IOFBF, sizeof(outBuf));

    printf("%d %d\n", (int)mesh.v.size(), (int)mesh.f.size());

    for (const Vec3& p : mesh.v) {
        printf("v %.15g %.15g %.15g\n", p.x, p.y, p.z);
    }

    for (const Tri& t : mesh.f) {
        printf("f %d %d %d\n", t.a + 1, t.b + 1, t.c + 1);
    }
}

int main() {
    startTime = chrono::steady_clock::now();

    readInput();

    bestMesh.v = origV;
    bestMesh.f = origF;

    renderResolution = (N > 180000 ? 160 : (N > 60000 ? 192 : 256));
    originalRender = renderMesh(origV, origF, renderResolution);

    if (N <= 40) {
        outputMesh(bestMesh);
        return 0;
    }

    tryEllipsoidCandidates();
    tryOrderedTorusGrid();
    tryOrderedSphereGrid();
    tryConservativeQEM();
    tryLocalVertexRemoval();

    outputMesh(bestMesh);
    return 0;
}
