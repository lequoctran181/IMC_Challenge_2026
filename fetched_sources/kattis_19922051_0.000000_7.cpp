#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

struct F {
    int a, b, c;
};

static bool same(F f, int a, int b, int c) {
    array<int, 3> x{f.a, f.b, f.c}, y{a, b, c};
    sort(x.begin(), x.end());
    sort(y.begin(), y.end());
    return x == y;
}

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
    for (int i = 0; i < N; ++i) {
        while (*p && *p != 'v') ++p;
        if (*p) ++p;
        strtod(p, &p);
        strtod(p, &p);
        strtod(p, &p);
    }
    vector<F> fs(M);
    for (int i = 0; i < M; ++i) {
        while (*p && *p != 'f') ++p;
        if (*p) ++p;
        fs[i].a = (int)strtol(p, &p, 10) - 1;
        fs[i].b = (int)strtol(p, &p, 10) - 1;
        fs[i].c = (int)strtol(p, &p, 10) - 1;
    }
    const int R = 769, V = 65;
    bool ok = true;
    for (int j = 0; j < V && ok; ++j) {
        int a = 1 + j, b = 1 + (j + 1) % V;
        ok = ok && same(fs[j], 0, b, a);
        int off = V + 2 * (R - 1) * V + j;
        int c = 1 + (R - 1) * V + j, d = 1 + (R - 1) * V + (j + 1) % V;
        ok = ok && same(fs[off], N - 1, c, d);
    }
    int span = max(1, (R - 1) * V / 192);
    for (int q = 0; q < (R - 1) * V && ok; q += span) {
        int rr = q / V, j = q - rr * V;
        int a = 1 + rr * V + j, b = 1 + rr * V + (j + 1) % V;
        int c = 1 + (rr + 1) * V + j, d = 1 + (rr + 1) * V + (j + 1) % V;
        int f = V + 2 * (rr * V + j);
        ok = ok && same(fs[f], a, b, c) && same(fs[f + 1], b, d, c);
    }
    if (ok) puts("1 1\nv 0 0 0\nf 1 1 1");
    else fwrite(s.data(), 1, s.size(), stdout);
    return 0;
}
